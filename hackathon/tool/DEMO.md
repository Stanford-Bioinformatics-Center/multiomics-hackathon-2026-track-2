# Demo script (under 4 minutes)

Setup before the talk (from `hackathon/tool/`, env `~/miniconda3/envs/mprobe` on danilogin):
- terminal 1: `python -m http.server --directory site 8000`, browser on http://localhost:8000
- terminal 2: `marimo run app.py` (store preloaded)
- terminal 3: empty, for the live commands

Numbers below are from the current build (`make all`); every one is in the reports' tables.

**0:00. Gallery, first row: all disease signatures side by side (`mprobe compare`).**
- Say: "Any disease signature in, MoTrPAC out, by layer and time. Blue = exercise pushes the disease genes back."
- The two mitochondria-down diseases are opposed in a median 54% of columns, the others in 4%:
  - type 2 diabetes (Mootha 2003) is opposed in 8 of 12 columns: human RNA EE 24 h t +10.7, rat gastrocnemius
    protein F 8 wk t +11.6;
  - PAH muscle (Malenfant 2015) in 5 of 12.
- Heart failure, aging (Liu 2013), septic shock and PAH blood are flat.
- Controls: Hostrup HIIT is red everywhere (same direction, t −9.1 in rat protein F 8 wk, as a positive control
  should be); the random abundance-matched set is flat.

**0:30. PAH muscle report.**
- Headline table:
  - human RNA at EE 24 h is opposed (15/18 genes, t +4.2);
  - human protein at RE 24 h goes the disease's way (2/19 opposed, t −4.1);
  - rat gastrocnemius protein is opposed by week 4 (F 4 wk 16/19, t +4.8).
- Section E sentence: "At EE–CON 24 h, RNA opposed 15/18 and protein 4/18; genome-wide RNA–protein agreement here is
  ρ = 0.09."
- Click a gene in the interactive grid: its trajectory in both layers.
- Section G: the caveat block and the guardrails table.

**1:15. The specificity lesson, live.**
- In the app, paste the ten-line random MitoCarta list from `examples/random_mito9.csv` (or run
  `mprobe run --signature examples/random_mito9.csv --name demo`, about 15 s).
- 7 of 9 random mitochondrial genes are opposed in rat protein F 8 wk: 99th–100th percentile of the abundance null,
  but only the 31st percentile of the pathway-class null.
- PAH at the same comparison: 91st percentile of the class null; type 2 diabetes: 98th.
- Say: "Exercise raises mitochondria. The class null tells you whether your disease adds anything beyond its
  pathway class."

**2:00. Discordance view, rat gastrocnemius** (gallery row "Layer discordance: rat SKM-GN", or
`mprobe discord --species rat --tissue SKM-GN`):
- timescale: RNA moves early, protein accumulates over weeks;
- agreement by timepoint;
- detection power: which layer responded is predicted with AUC 0.67 from detection alone (baseline mRNA,
  proteomics missing values) vs 0.73 for the full model (shuffled labels 0.49). Say: "much of 'which layer
  responded' is which assay could see the gene."
- the three-layer pathway panel: TCA and other pathways in RNA, protein and metabolites (pools, not flux).

**3:00. Library.**
- `mprobe library query --signature examples/pah_muscle_malenfant2015.csv` ranks every MoTrPAC contrast:
  - top opposed: rat vastus lateralis RNA, female 4 wk;
  - top concordant: rat brown adipose RNA, male 8 wk.
- `store/library/motrpac_library.gmt` loads into Enrichr-style tools.

**3:30. README guardrails table.** Safe statement vs unsafe upgrade. End.
