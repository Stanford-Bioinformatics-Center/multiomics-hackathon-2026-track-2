"""Save every explorer example's result as static JSON for the website.

The website loads these files directly when an example is chosen, so the examples work as a
backup even if the analysis API is not running. Uploaded lists always run the live analysis.
Run from the repository root after changing the explorer, the examples or the MoTrPAC data:

    hackathon/tool/.venv/Scripts/python.exe apps/api/scripts_save_examples.py
"""
import datetime as dt
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from motrpac_probe_service import explorer  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "web" / "public" / "examples"
OUT.mkdir(parents=True, exist_ok=True)
stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
index = []
for item in explorer.list_examples():
    result = explorer.analyse([{"kind": item["kind"], "example": item["id"]}])
    result = {**result, "saved": {"generated_at": stamp, "example": item["id"],
                                  "note": "Saved result of the live explorer analysis for this example."}}
    text = json.dumps(result, allow_nan=False, separators=(",", ":"))
    (OUT / f"{item['id']}.json").write_text(text, encoding="utf-8")
    index.append({**item, "file": f"{item['id']}.json", "generated_at": stamp})
    print(f"{item['id']}: {len(text) / 1024:.0f} KB")
(OUT / "index.json").write_text(json.dumps({"examples": index}, indent=1), encoding="utf-8")
print("saved to", OUT)
