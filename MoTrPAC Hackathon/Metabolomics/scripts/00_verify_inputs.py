"""Check the bundled input files against data/source_manifest.json."""

import hashlib
import json

from common import DATA

manifest = json.loads((DATA / "source_manifest.json").read_text(encoding="utf-8"))
for name, entry in manifest.items():
    path = DATA / entry["file"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise ValueError(f"SHA-256 mismatch for {entry['file']}")
    print(f"Verified {entry['file']}")

factors = json.loads((DATA / manifest["st000763_factors"]["file"]).read_text(encoding="utf-8"))
data = json.loads((DATA / manifest["st000763_data"]["file"]).read_text(encoding="utf-8"))
if len(factors) != 218 or len(data) != 673:
    raise ValueError(f"Expected 218 samples and 673 features; found {len(factors)} and {len(data)}")
print("Verified 218 ST000763 samples and 673 metabolite features.")
