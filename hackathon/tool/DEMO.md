# Demo script (under 4 minutes)

Setup before the talk (from `hackathon/tool/`, env `~/miniconda3/envs/mprobe` on danilogin):
- terminal 1: `python -m http.server --directory site 8000`, browser on http://localhost:8000
- terminal 2: `marimo run app.py` (store preloaded)
- terminal 3: empty, for the live commands

Numbers below are from the current build (`make all`); every one is in the reports' tables.

**0:00. The app, the team's PAH analysis** (`marimo run app.py`; default example "PAH, 9 lower muscle proteins").
- Time course tab: the team's rat plot, generalised. All nine proteins rise above sedentary by 8 weeks of training
  (6 of 9 at FDR < 0.05, orange), while RNA (blue) barely moves for most: protein and RNA disagree.
- "Do RNA and protein agree?": per gene, protein-only responses; the detection columns show every gene is well
  measured in both assays, so this is not a detection artefact.
- Switch "Exercise data" to Human: after one bout it is the reverse, RNA up at 24 h and protein flat.
- Drag the FDR slider to 0.10: stars and classes update live.

**1:15. Generalisable: pick another disease** in the Example menu, e.g. "Type 2 diabetes, skeletal muscle
(Mootha 2003)": same page, 87 genes, 12 panels; or paste any list (`GENE,-1` per line).
- "Is it specific?": each cell shows cameraPR t and the class-null percentile. PAH's nine proteins sit around the
  75th–95th percentile of random mitochondrial sets: exercise raises mitochondria in general.

**2:00. Several diseases at once** (gallery first row, `mprobe compare`): mitochondria-down signatures (PAH, type 2
diabetes) are opposed, heart failure / aging / sepsis are not; the positive control (HIIT proteome) moves the same way.

**2:40. Discordance report** (gallery row "Layer discordance: rat SKM-GN", `mprobe discord`): per timepoint the
share of each layer that changes, RNA–protein ρ, and the detection-power result (AUC 0.67 from detection alone vs
0.73 full).

**3:20. Library and the team's blood analysis**: `mprobe library query ...` ranks every MoTrPAC contrast;
`mprobe run --signature examples/pah_blood_gse33463_ranked.csv.gz` reproduces the team's blood rank correlations and
27 / 13 GO:BP concordance exactly. Close on the README guardrails table.
