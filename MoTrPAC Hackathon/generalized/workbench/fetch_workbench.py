"""Fetch one named Metabolomics Workbench study by validated ST accession.

Only the official /factors and named-metabolite /data REST endpoints are used.
The separate workbench_adapter.py requires a reviewed study design before any
statistical contrast is calculated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen


BASE = "https://www.metabolomicsworkbench.org/rest/study/study_id"
MAX_BYTES = 100 * 1024 * 1024


def validate_accession(value: str) -> str:
    accession = value.upper()
    if not re.fullmatch(r"ST\d{6}", accession):
        raise ValueError("Study accession must have the form ST followed by six digits")
    return accession


def fetch_json(accession: str, endpoint: str) -> tuple[bytes, str]:
    url = f"{BASE}/{accession}/{endpoint}"
    with urlopen(url, timeout=60) as response:
        final_url = response.geturl()
        parsed = urlparse(final_url)
        if parsed.scheme != "https" or parsed.hostname not in {
            "www.metabolomicsworkbench.org", "metabolomicsworkbench.org"
        }:
            raise ValueError(f"Unexpected redirect away from Metabolomics Workbench: {final_url}")
        payload = response.read(MAX_BYTES + 1)
    if len(payload) > MAX_BYTES:
        raise ValueError(f"{endpoint} response exceeds {MAX_BYTES} bytes")
    parsed_json = json.loads(payload)
    records = list(parsed_json.values()) if isinstance(parsed_json, dict) else parsed_json
    if not isinstance(records, list) or not records or not all(isinstance(x, dict) for x in records):
        raise ValueError(f"{endpoint} endpoint did not return a nonempty record collection")
    if any(x.get("study_id") != accession for x in records):
        raise ValueError(f"{endpoint} endpoint returned records from another study")
    return payload, url


def download_study(accession: str, output_root: Path, overwrite: bool = False) -> Path:
    accession = validate_accession(accession)
    target = output_root / accession
    paths = {name: target / f"{name}.json" for name in ("factors", "data")}
    manifest_path = target / "source_manifest.json"
    if not overwrite and any(path.exists() for path in [*paths.values(), manifest_path]):
        raise FileExistsError(f"{target} already has Workbench files; pass --overwrite to replace them")
    # Fetch and validate both payloads before writing either file.
    downloaded = {name: fetch_json(accession, name) for name in paths}
    target.mkdir(parents=True, exist_ok=True)
    manifest = {
        "study_id": accession,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "source": "Metabolomics Workbench REST API, named-metabolite study data",
        "files": {},
    }
    for name, path in paths.items():
        payload, url = downloaded[name]
        path.write_bytes(payload)
        manifest["files"][name] = {
            "path": path.name, "url": url, "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accession", help="Workbench study accession, e.g. ST000763")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Parent directory for the accession subfolder")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing downloaded files")
    args = parser.parse_args()
    target = download_study(args.accession, args.output_dir, args.overwrite)
    print(f"Saved named Workbench factors/data and SHA256 manifest in {target}")


if __name__ == "__main__":
    main()
