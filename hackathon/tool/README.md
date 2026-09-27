# motrpac_probe: how much of a disease signature does exercise oppose, by omic layer and time?

`mprobe` takes any disease signature (genes or proteins up/down in disease; metabolites for the METAB layer) and
reports, for every MoTrPAC exercise comparison, whether the exercise response moves those molecules against the
disease direction ("opposed") or with it, separately for RNA, protein, phosphosites and metabolites, at every
timepoint. It also reports why the answer differs between layers, and how much of the answer is the signature's
pathway class rather than the disease. It uses only public MoTrPAC summary statistics.

## Read these first

| file | what it shows |
|---|---|
| `analysis.py` | the whole story in ~60 commented lines: load MoTrPAC, score a signature by layer and time, check specificity, compare RNA with protein (`python analysis.py`, 2 s) |
| `src/motrpac_probe/core.py` | the scoring: does exercise move each gene against the disease, and is the whole set shifted (cameraPR) |
| `src/motrpac_probe/nulls.py` | specificity: the same test on random sets of the same abundance or pathway class |
| `src/motrpac_probe/layers.py` | which layer answered, next to whether it could have (detection power) |
| `src/motrpac_probe/ranked.py` | the team's blood analysis, generalised: whole-ranking correlation and pathway concordance |
| `src/motrpac_probe/store.py` | how the MoTrPAC tables become one gene-level store |

Everything else is plumbing: report and figure rendering (`run.py`, `figures.py`, `render.py`), the CLI, the gallery
and the app. Method details for every module: `docs/METHODS.md`.

## Quick start (macOS, Linux, Windows)

Paste the blocks as they are (they contain no comments, so zsh and PowerShell accept them). Nothing needs R, a
compiler or cluster access.

**1. Clone and create the environment** (recommended: conda or mamba; the same commands on every OS).

```bash
git clone -b t3code/build-motrpac-probe-tool https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-2.git
cd multiomics-hackathon-2026-track-2/hackathon/tool
conda env create -f environment.yml
conda activate mprobe
```

Without conda: a plain virtual environment with Python 3.11 to 3.14 works too (on Windows PowerShell use
`.venv\Scripts\Activate.ps1` instead of the `source` line).

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt -e .
```

**2. Data.** Nothing to do: the first command that needs data downloads the store once (about 180 MB of public
MoTrPAC summary statistics from the GitHub release
[`mprobe-store-v1`](https://github.com/Stanford-Bioinformatics-Center/multiomics-hackathon-2026-track-2/releases/tag/mprobe-store-v1),
every file sha256-checked). To do it explicitly, or from a local copy: `mprobe store fetch [--url <path-or-URL>]`.

**3. See the core result in the terminal** (about 2 s): genes opposed / measured, cameraPR t, the specificity
percentile and RNA vs protein, for the PAH muscle signature and then for type 2 diabetes.

```bash
python analysis.py
python analysis.py examples/type2_diabetes_muscle_mootha2003.csv
```

**4. Full reports from the command line.** Each command writes `out/<name>/report.html`; open it in a browser
(`open out/pah_muscle_malenfant2015/report.html` on macOS). In order: a gene / protein signature (~15 s); random
mitochondrial genes (the specificity lesson); a full disease ranking (sections R and P, the team's blood analysis);
a metabolite signature; several diseases side by side; which omic layer responds and why (~1 min); MoTrPAC as GMT
gene sets and a query against them; and a gallery of every report at http://localhost:8000.

```bash
mprobe run --signature examples/pah_muscle_malenfant2015.csv
mprobe run --signature examples/random_mito9.csv --name demo
mprobe run --signature examples/pah_blood_gse33463_ranked.csv.gz
mprobe run --signature examples/tca_intermediates_demo.csv
mprobe compare --signatures examples/pah_muscle_malenfant2015.csv examples/type2_diabetes_muscle_mootha2003.csv examples/random_mito9.csv
mprobe discord --species rat --tissue SKM-GN
mprobe library build
mprobe library query --signature examples/pah_muscle_malenfant2015.csv
mprobe site build
python -m http.server --directory site 8000
```
Every option has a flag (`mprobe run -h`); `config.example.yaml` mirrors them all (`--config file.yaml`).

**5. Your own signature.** A CSV with one row per gene and a direction (+1 up in disease, −1 down), e.g.
`my_disease.csv`:

```
gene_symbol,direction
NDUFA9,-1
LDHA,1
```

```bash
mprobe run --signature my_disease.csv --name my_disease
```

**6. Interactive app (optional).** Everything above works without it. It needs one extra package; it opens
http://localhost:2718 and the store loads in about 4 s.

```bash
pip install "marimo>=0.16"
marimo run app.py
```
1. "Signature from": an example, an uploaded CSV, or pasted lines such as `NDUFA9,-1`.
2. Choose tissues, FDR cutoff, background and the number of random sets; press **Run analysis** (~8 s).
3. Read down the page: coverage and unmapped genes; the agreement grid with a gene selector for its trajectory;
   set-level opposition with null percentiles; RNA vs protein; the three-layer pathway panel; caveats.
4. **Export report** writes the same HTML report as `mprobe run`.

**7. Tests** (about 3 minutes; downloads the store first if needed). On a fresh clone, tests that need the
maintainers' raw hackathon files (golden numbers from the original analysis) skip automatically; everything else
runs against the downloaded store. No failures expected.

```bash
python -m pytest -q
```

**Maintainers** with the raw hackathon data (`hackathon/data/`, the MoTrPAC R packages) rebuild the store with
`mprobe store build`, everything else with `python scripts/build_all.py` (or `make all`; `make slurm-all` on the
cluster), and a new bundle with `mprobe store bundle`.

### Input formats

**Signature CSV:** `gene_symbol` (human; rat symbols, UniProt and Ensembl ids are also accepted, in the same column
or in `uniprot` / `ensembl` / `rat_symbol` columns), `direction` (+1 up in disease, −1 down; `up`/`down` accepted),
optional `group`, `weight`, `source`. Rat genes map to human through the RGD ortholog table used by
`scripts/03_join.py`. Every input row is reported with its mapping status (mapped, not found, conflicting
direction, duplicate). Rows in the group `phenotype` (configurable: `--context-groups`) are shown but not counted.
Weights are carried and displayed; the statistics are unweighted.

**Full disease ranking** (a genome-wide disease-vs-control contrast): a CSV with a gene column and a `t` (or
`stat`, `z`, `score`, `logFC`) column and no `direction` column, e.g. a limma table. `mprobe run` detects it and
adds sections R and P below; a directional signature derived from it (BH < 0.05, at most 250 genes per direction)
drives the other sections. `--exact-symbols` matches gene symbols exactly, as the team pipeline does.

**Metabolite signature:** `refmet_name` (or `hmdb` / `kegg`) and `direction`; scored on the METAB layer.

## Architecture (the team's whiteboard design)

| whiteboard box | module |
|---|---|
| Frontend: dynamic interactive charts, toggles, overlay | reports (inline SVG + Vega-Lite grid with tissue / layer filters and click-a-gene trajectories), `app.py` (marimo), `site/` gallery |
| User input → Query (gene / protein set) → parser | `signature.py` (symbols, UniProt, Ensembl, rat symbols, pasted lists), `ranked.py` (full disease rankings), `metab.py` (RefMet / HMDB / KEGG) |
| MoTrPAC → transcriptomics / proteomics / metabolomics DB → table join | `store.py`: one parquet store, 422 gene-level columns + 217 metabolite columns, joined on human gene symbol / RefMet name |
| Flag | detection-power flags, human 15–45 min flag, fibre-type markers, caveat block (`layers.py`, `caveats.py`) |
| Key/value tags, processing, plots per layer | `provenance.json` + parameters; `core.py`, `nulls.py`, `everywhere.py`, `discord.py`; `figures.py`, `render.py` |
| Next steps: other omics, discordance associated with diseases | `compare` (several diseases side by side), `discord` (which layer responds and why), `library` (GMT for any tool) |
| Optional LLM | not included: every report sentence is filled in from computed numbers, so results are reproducible |

## Relationship to the team pipeline on main (`MoTrPAC Hackathon/`)

The team's PAH pipeline (15 numbered R/Python scripts, collection c2.0) is one case study; this tool is its
generalisation to any disease, tissue, species and layer, and reproduces it (`tests/test_team.py`, reference values
are the team's committed outputs in `tests/fixtures/team_main_c2/`). Both use MotrpacHumanPreSuspensionAnalysis
2.0.8, collection c2.0; the exported tables are identical row by row.

| team script | what it does | in the tool | reproduced |
|---|---|---|---|
| 00, 03–04 | nine lower PAH proteins: UniProt match; RNA at three times | `run` sections A–C (any signature; UniProt, Ensembl and rat ids accepted) | 0/9, 7/9, 9/9 positive; 7 at BH < 0.05 at 24 h; 0/27 protein |
| 04.5–05 | OXPHOS membership; MoTrPAC's precomputed OXPHOS test | cameraPR port validated against the package's CAMERA results; OXPHOS/TCA pool null | package values to 1e-8 |
| 06–07 | GSE33463 IPAH-minus-healthy PBMC limma ranking | stays in the team pipeline; its output is the example `pah_blood_gse33463_ranked.csv.gz` | input |
| 07.5, 08 | MoTrPAC blood gene ranks; Spearman of disease vs exercise ranking, gene-label permutation | section R: every comparison (422), Fisher 95% CI, BH, calibration against non-exercise reference contrasts | ρ to 3 decimals at all six times |
| 09–11 | GO:BP cameraPR on the disease ranking; join to MoTrPAC's precomputed GO:BP results | section P: same test on both sides for every comparison, plus the team's precomputed-results variant | 5,183 sets, 63 at BH < 0.05; 40 joint rows = 27 same + 13 opposite |
| 12 | blood sensitivity checks | the sensitivity section (FDR, background, null size) | — |

Conventions: ρ > 0 means exercise moves genes the same way as the disease (as in the team pipeline); cameraPR
t > 0 in sections B–D means exercise opposes a directional signature. Both are labelled in every table.

## What a report contains (`mprobe run`)

The first screen is the headline table: for each tissue × layer × time, genes opposed / measured, the signed
cameraPR t, and the pathway-class null percentile, then one coverage line and one caveat line. Sections follow:

| | section | method |
|---|---|---|
| A | coverage and detection power | which genes are measured where; RNA baseline tertile, proteomics missing values, GTEx v8 TPM in the matched tissue, features collapsed |
| B | agreement grid | cell = sign(exercise logFC) × disease direction × min(\|stat\|, 4); −1 = exercise opposes; dot = fdr_bh < 0.05; per column n measured, n opposed, binomial sign test, and the same test at the effective number of genes (Nyholt; NARRATIVE R3) |
| C | set-level opposition (primary) | cameraPR (validated Python port of limma) on the full ranking with member statistics multiplied by −direction, so t > 0 = opposed; up- and down-halves separately; BH FDR |
| D | calibration nulls | (i) 1,000 random sets matched gene-by-gene on abundance decile; (ii) 1,000 sets matched on MitoCarta / GO:CC complex / secreted class ("your signature or its pathway class?"); optional fixed-pool null (NARRATIVE R1a) |
| E | layer discordance | for genes measured in RNA and protein at the same tissue × time: opposed in each layer, signature vs genome-wide Spearman ρ, "which layer answered" per gene next to its detection flags, one templated sentence per comparison |
| F | everywhere | the same scores in all 422 MoTrPAC comparison columns (19 rat tissues, human muscle/adipose/blood, incl. non-exercise reference contrasts as a negative control), ranked by opposition and by same-direction |
| G | caveats | always printed: human 15–45 min flag, fibre-type markers, optimistic sign tests, cross-species/tissue mismatch, and the guardrails table below |
| | sensitivity | every conclusion re-derived at FDR 0.10, muscle-intrinsic background, 200 vs 1,000 null sets; flips are listed |

All prose in reports is filled in from computed numbers; there is no free-text generation. Every number in prose is
also in a table (`tables/*.csv`); figures are PNG at native aspect ratio plus inline SVG. `provenance.json` records
the git SHA, package versions, store hash, signature hash, parameters, timings and the command.

Definitions are the ones in `../NARRATIVE.md`. The code that computes them is executed from the analysis scripts
themselves (`legacy.py` lifts `sign_test`, `camera_pr`, `bh`, `collapse`, `parse_human`, `rho_ci`, `boot_rho` …
from `scripts/03, 05, 06, 08, 10` without running those pipelines), and the tests reproduce the deck numbers.

## Interpretation guardrails

From the team's planning workbook (`MoTrPAC_PAH_Comparison_Analysis_Map.xlsx`, sheet "Start Here", section 4),
generalised to any disease signature, plus two rows from our own results. Printed in every report.

| Safe statement | Unsafe upgrade |
|---|---|
| A disease-signature pathway is directionally opposed after exercise. | Exercise reverses or treats the disease. |
| Blood and muscle share a pathway label. | The same genes changed or one tissue signalled to the other. |
| An early phosphosite and later target-motif signal are temporally compatible. | The phosphosite caused the later transcriptional response. |
| A plasma protein matches tissue RNA/protein in time and direction. | The tissue is proven to be the source of the plasma protein. |
| A PLIER LV is associated with a pathway prior. | The pathway caused the LV or independently validates it. |
| RNA for a gene rose after one bout. | The protein rose, or will rise. |
| A gene is "protein-only" in one contrast. | The gene is post-transcriptionally regulated. |

Further rules the reports enforce: the sign test is descriptive (genes are co-regulated), cameraPR is the evidence;
read the class-null percentile before attributing opposition to the disease; human 15–45 min columns partly measure
biopsy composition (NARRATIVE R2); contractile signatures must be read against fibre-type markers (rat training
shifts IIb → IIa/IIx).

## Why this exists

What already exists: the MoTrPAC Data Hub and its visualisations; the R packages MotrpacRatTraining6moData and
MotrpacHumanPreSuspensionAnalysis (all differential-analysis tables, with enrichment helpers); the Enrichr library
`MoTrPAC_2023` (225 gene sets, https://maayanlab.cloud/Enrichr/). This tool's own GMT covers all 422 human and rat
comparison columns, with direction and provenance. What this adds: (1) layer × time resolution for a user's own signature, RNA and protein
side by side; (2) calibration nulls including a pathway-class null, so "opposed" can be checked against "any
mitochondrial list looks opposed"; (3) detection-power flags next to every "which layer responded" call; (4) a
metabolite layer at pathway level; (5) a fixed caveat block and positive/negative/pathway-class controls; (6)
MoTrPAC as a GMT library with provenance. What it does not do: per-person or per-animal analysis (summary
statistics only), causal or therapeutic inference, flux, cell-type deconvolution, or any re-analysis of raw data.

## Customizing

- Colours and fonts: `src/motrpac_probe/templates/style.css` (HTML) and the constants at the top of `render.py`
  (figures). Report sections can be switched off with `--sections` / `sections:` in the config.
- New gene-set collections: add rows (gs_name, gene_symbol, source) to `store/genesets.parquet` or rebuild with an
  extra CSV in `store.build_genesets`; any set is usable for the fixed-pool null (`--pool-sets`).
- New pathway maps for the metabolite panel: edit `store/pathway_map.csv` (pathway, side, member, rule, source).
- New examples: add a CSV to `examples/` and an entry to `examples/examples.json`; `make all` picks it up.

## Signature library

`mprobe library build` writes every MoTrPAC comparison (422 columns) as UP and DOWN gene sets (fdr_bh < 0.05,
ranked by |stat|, at most 250 genes) to `store/library/motrpac_library.gmt`, with `motrpac_library_index.json`
(tissue, layer, species, sex, time, n, provenance). GMT is the Enrichr / GSEA library format, so the file loads into
Enrichr-style tools. `library query` ranks all columns by set-level opposition and concordance with a signature,
with Jaccard overlaps.

## Data sources

Only public summary statistics are used; no controlled-access or participant-level data, so no data-use agreement
beyond the public MoTrPAC terms applies.
- MoTrPAC Data Hub: https://motrpac-data.org
- MotrpacRatTraining6moData 2.0.0 (rat endurance training, 19 tissues): https://github.com/MoTrPAC/MotrpacRatTraining6moData
- MotrpacHumanPreSuspensionAnalysis 2.0.8 (human acute exercise, pre-suspension release): https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis
- Rat metabolomics also deposited at Metabolomics Workbench, project PR001020; RefMet: https://www.metabolomicsworkbench.org
- GTEx v8 median gene TPM (CFDE): https://gtexportal.org
- MitoCarta3.0: Rath et al. 2021 Nucleic Acids Res; GO:CC and pathway sets via msigdbr (MSigDB)
- Malenfant et al. 2015 J Mol Med (PAH vastus lateralis proteome); Cheadle et al. 2012 (GSE33463, PAH blood
  erythroid signature); Hostrup et al. 2022 eLife 11:e69802 (HIIT muscle proteome). Other example sources are cited
  in `examples/README.md`.

## Limitations

Group-level contrasts only; human protein and phospho have no per-feature missingness in the public tables; rat
samples are ~48 h after the last bout; the human acute study has one bout and no training time course; signature
quality (often 4 vs 4 patients) limits every conclusion; cameraPR assumes an inter-gene correlation of 0.01 and the
nulls treat genes as exchangeable within strata; the pathway-class null has only three binary class flags.

## Tests

`make test` (pytest): golden numbers from `../deck_extracts/deck_numbers.csv` (15/18, 2/19, 15/19, down-half t,
OXPHOS full vs muscle-intrinsic background, specificity percentiles 97.7/82/77.8, R1b, effective n, RNA–protein ρ,
the 258-unit everywhere summary), cameraPR vs the package's precomputed CAMERA results, the store schema, signature
mapping, CLI smoke runs on every example with report validation, and `tests/test_matrix.py` (robustness matrix;
summary in `tests/REPORT.md`).
