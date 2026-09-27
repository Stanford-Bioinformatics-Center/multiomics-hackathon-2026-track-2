# Regenerating the report gallery from a clean commit

The gallery reports committed under `hackathon/tool/site/` were generated on an earlier commit
(`5e2cde3`) that is not on the `integration/exercise-signature-explorer` branch. Their
`provenance.json` therefore records a stale `repo_git_sha`. The generation pipeline itself is
correct: regenerating on this branch records the current commit (verified — a fresh
`mprobe run` records the integration HEAD, e.g. `a1dbfc4`, not `5e2cde3`).

To regenerate the full gallery from a clean checkout (do this on a clean working tree so provenance
does not record a dirty state, or accept the recorded dirty flag):

```sh
cd hackathon/tool
. .venv/bin/activate            # or a fresh: python -m venv .venv && pip install -e ".[test,app]"
mprobe store fetch              # if the store is not present
python scripts/build_all.py     # or: make all   (regenerates every example report + the gallery)
mprobe site build               # rebuild site/index.html
```

Then confirm the recorded SHA matches HEAD:

```sh
python - <<'PY'
import json, subprocess
head = subprocess.run(["git","rev-parse","HEAD"],capture_output=True,text=True).stdout.strip()
p = json.load(open("site/pah_muscle_lower9_malenfant2015/provenance.json"))
assert p["repo_git_sha"] == head, (p["repo_git_sha"], head)
print("gallery provenance SHA matches HEAD:", head)
PY
```

Note: regenerating overwrites teammate-authored committed reports, so it should be done as its own
reviewed commit, not mixed with feature work.
