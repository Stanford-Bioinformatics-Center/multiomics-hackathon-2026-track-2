"""Run one or more processed-signature queries from an app-friendly JSON file.

The runner is deliberately small: source-specific preprocessing belongs in an
adapter, while query_core and plots own the scientific operations. Each job is
independent so an unsupported context can be reported without hiding the
results of compatible contexts.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

from query_core.engine import run_query
from query_core.schema import SPECIES


JOB_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
FILTERS = {
    "species", "reference_species", "tissue", "disease_contrast",
    "reference_contrast_category", "reference_contrast", "timepoints",
    "rank_field", "q_threshold", "min_abs_log2_fc", "ortholog_map",
}


def _source_path(config_dir: Path, value: str, field: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty path")
    path = Path(value)
    path = path if path.is_absolute() else config_dir / path
    # Source availability is checked while running each job. One missing
    # upload must not prevent other valid jobs in the same config from running.
    return path.resolve()


def load_jobs(config_path: Path) -> list[dict]:
    config_path = Path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise ValueError("Config must be an object with schema_version: 1")
    jobs = config.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("Config needs a nonempty jobs array")
    seen: set[str] = set()
    resolved = []
    for i, job in enumerate(jobs, 1):
        if not isinstance(job, dict):
            raise ValueError(f"jobs[{i}] must be an object")
        job_id = job.get("id")
        if not isinstance(job_id, str) or not JOB_ID.fullmatch(job_id):
            raise ValueError(f"jobs[{i}].id must be a safe 1–64 character ID")
        if job_id in seen:
            raise ValueError(f"Repeated job ID: {job_id}")
        seen.add(job_id)
        filters = job.get("filters", {})
        if not isinstance(filters, dict) or set(filters) - FILTERS:
            raise ValueError(f"jobs[{i}].filters has unsupported fields: {sorted(set(filters) - FILTERS) if isinstance(filters, dict) else 'non-object'}")
        filters = dict(filters)
        if "ortholog_map" in filters:
            filters["ortholog_map"] = _source_path(config_path.parent,
                                                    filters["ortholog_map"],
                                                    f"jobs[{i}].filters.ortholog_map")
        if "timepoints" in filters:
            times = filters["timepoints"]
            if (not isinstance(times, list) or
                    not all(isinstance(x, str) and x.strip() for x in times)):
                raise ValueError(f"jobs[{i}].filters.timepoints must be an array of nonempty strings")
            filters["timepoints"] = tuple(times)
        for field in ("species", "reference_species", "tissue", "disease_contrast",
                      "reference_contrast_category", "reference_contrast"):
            if field in filters and (not isinstance(filters[field], str) or
                                     not filters[field].strip()):
                raise ValueError(f"jobs[{i}].filters.{field} must be a nonempty string")
        for field in ("species", "reference_species"):
            if field in filters and filters[field].lower() not in SPECIES:
                raise ValueError(f"jobs[{i}].filters.{field} must be human or rat")
        if "rank_field" in filters and filters["rank_field"] not in {
                "auto", "log2_fc", "statistic"}:
            raise ValueError(f"jobs[{i}].filters.rank_field is unsupported")
        for field in ("q_threshold", "min_abs_log2_fc"):
            if field in filters:
                value = filters[field]
                if (isinstance(value, bool) or not isinstance(value, (int, float)) or
                        not math.isfinite(value) or
                        (field == "q_threshold" and not 0 < value <= 1) or
                        (field == "min_abs_log2_fc" and value < 0)):
                    raise ValueError(f"jobs[{i}].filters.{field} is outside its valid range")
        plots = job.get("plots", False)
        if isinstance(plots, bool):
            plots = {"enabled": plots, "top": 25}
        elif isinstance(plots, dict):
            if set(plots) - {"enabled", "top"}:
                raise ValueError(f"jobs[{i}].plots has unsupported fields")
            plots = {"enabled": plots.get("enabled", True), "top": plots.get("top", 25)}
        else:
            raise ValueError(f"jobs[{i}].plots must be Boolean or object")
        if not isinstance(plots["enabled"], bool) or type(plots["top"]) is not int or plots["top"] < 1:
            raise ValueError(f"jobs[{i}].plots needs Boolean enabled and positive integer top")
        resolved.append({
            "id": job_id,
            "disease": _source_path(config_path.parent, job.get("disease"), f"jobs[{i}].disease"),
            "reference": _source_path(config_path.parent, job.get("reference"), f"jobs[{i}].reference"),
            "filters": filters,
            "plots": plots,
        })
    return resolved


def run_config(config_path: Path, out_dir: Path) -> dict:
    config_path = Path(config_path).resolve()
    jobs = load_jobs(config_path)
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "config_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(),
        "jobs": [],
    }
    for job in jobs:
        job_dir = out_dir / job["id"]
        record = {
            "id": job["id"], "status": "error",
            "disease": str(job["disease"]), "reference": str(job["reference"]),
            "results_dir": str(job_dir), "filters": {
                k: list(v) if isinstance(v, tuple) else str(v) if isinstance(v, Path) else v
                for k, v in job["filters"].items()
            },
        }
        try:
            summary = run_query(job["disease"], job["reference"], job_dir,
                                **job["filters"])
            record["status"] = "complete"
            record["counts"] = summary["counts"]
            record["files"] = {name: str(job_dir / name) for name in (
                "summary.json", "matches.csv.gz", "coverage.csv", "unmatched.csv",
                "ambiguous.csv.gz", "summary_by_group.csv", "rank_correlation.csv",
            )}
            if job["plots"]["enabled"]:
                from plots.plot_results import render
                plot_dir = job_dir / "plots"
                render(job_dir, plot_dir, top=job["plots"]["top"])
                record["plots_dir"] = str(plot_dir)
        except Exception as exc:
            # Jobs are independent. Preserve the failure for the app while
            # allowing other requested contexts to finish.
            record["status"] = "error"
            record["error"] = f"{type(exc).__name__}: {exc}"
        manifest["jobs"].append(record)
        pending = out_dir / "manifest.json.tmp"
        pending.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        pending.replace(out_dir / "manifest.json")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = run_config(args.config, args.out)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        parser.error(str(exc))
    n_ok = sum(job["status"] == "complete" for job in manifest["jobs"])
    print(f"Completed {n_ok}/{len(manifest['jobs'])} jobs; manifest: {Path(args.out).resolve() / 'manifest.json'}")
    return 0 if n_ok == len(manifest["jobs"]) else 1


if __name__ == "__main__":
    sys.exit(main())
