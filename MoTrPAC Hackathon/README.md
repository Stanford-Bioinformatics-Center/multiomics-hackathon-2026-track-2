# MoTrPAC–PAH hackathon analysis

This repository compares **resting PAH measurements** with **acute endurance-exercise results in a separate, generally healthy MoTrPAC cohort**. It has two modules: skeletal-muscle proteins/RNA and PBMC/whole-blood RNA. The comparisons generate hypotheses; they do not test whether exercise treats PAH.

## Run it

Use Python 3.10 or newer and R 4.4 or newer. The Python scripts use only the standard library. The R scripts require `limma`, `dplyr`, and `MotrpacHumanPreSuspensionAnalysis` **0.2.4**. The [MoTrPAC package installation guide](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis#installation) explains which Bioconductor version goes with each R version. The local verification used R 4.4.2, limma 3.62.2, and package 0.2.4. The local package source checkout is MoTrPAC commit `17a60dbe2e506c8ae912653a7bf979ff35656e6d` (20 July 2026). After setting up Bioconductor as described in that guide, the same source can be installed in R with:

```r
BiocManager::install("limma")
pak::pak("MoTrPAC/MotrpacHumanPreSuspensionAnalysis@17a60dbe2e506c8ae912653a7bf979ff35656e6d")
```

Check the installed version with:

```sh
Rscript -e 'stopifnot(as.character(packageVersion("MotrpacHumanPreSuspensionAnalysis")) == "0.2.4")'
```

From this folder:

```sh
python scripts/run_pipeline.py
python scripts/run_pipeline.py --module muscle
python scripts/run_pipeline.py --module blood
```

The runner finds `Rscript` on `PATH` or in a standard Windows R installation. Use `--rscript /path/to/Rscript` if needed. Each numbered script also runs by itself. Run dependencies in the order below. The runner stops on the first error; it never downloads data or installs packages.

The required GEO matrix, platform annotation, and curated PAH protein table are already in `data/`. The package supplies published MoTrPAC summaries and GO definitions. Script 07.5 regenerates `data/processed/blood_gene_ranked.csv.gz`; it is an input to scripts 08 and 12. The existing result files are grouped under `outputs/PAH protein/` and `outputs/PAH blood/` and are overwritten by a rerun.

| Script | What it does | Main product |
| --- | --- | --- |
| 00 | Audits the nine matching muscle RNA estimates at **all three** sampled times | `outputs/PAH protein/pah_muscle_rna_timepoint_summary.csv` |
| 01–02 | Export published MoTrPAC muscle RNA and protein summary results; no model fit | `data/processed/motrpac_muscle_*.csv.gz` |
| 03–04 | Join nine PAH paper proteins to MoTrPAC protein by UniProt, then to 24-hour RNA by current gene symbol | `outputs/PAH protein/pah_*matches.csv` |
| 04.5–05 | Audit nine-gene GO OXPHOS membership and export MoTrPAC's precomputed broad pathway test | `outputs/PAH protein/*oxphos*.csv` |
| 06 | Read the submitter-processed GSE33463 PBMC matrix; select IPAH and healthy samples, check values | `outputs/PAH blood/gse33463_prepared.rds` |
| 07 | Fit an IPAH-minus-healthy **probe-level** limma model and map probes to symbols | `outputs/PAH blood/pah_blood_{probe,gene}_ranked.csv.gz` |
| 07.5 | Rebuild the published MoTrPAC blood gene ranks with the source project's feature-collapse rule | `data/processed/blood_gene_ranked.csv.gz` |
| 08 | Compare PAH and exercise gene ranks at six blood time points | `outputs/PAH blood/pah_motrpac_blood_rank_association_by_time.csv` |
| 09 | Test GO Biological Process sets against the PAH PBMC ranking | `outputs/PAH blood/pah_gobp_camera.csv` |
| 10–11 | Export MoTrPAC's blood pathway tests and join exact GO names to PAH results | `outputs/PAH blood/pah_motrpac_blood_gobp_*.csv*` |
| 12 | Check the blood rank correlation after two specific analytic changes | `outputs/PAH blood/pah_motrpac_blood_rank_sensitivity.csv` |

**Input provenance.** `data/from paper/pah_lower_proteins_malenfant2015.csv` is a transcription of Table 2's nine lower proteins from [Malenfant et al. (2015)](https://doi.org/10.1007/s00109-014-1244-0). Its `pah_to_control_ratio` is measured protein abundance in resting PAH muscle divided by control abundance: 0.71 for NDUFA9 corresponds to about 29% lower measured abundance. The paper abstract reports “9 downregulated proteins in PAH skeletal muscles.” `data/processed/GSE33463/GSE33463_series_matrix.txt.gz` and `data/processed/GPL6947.annot.gz` are the [GEO GSE33463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33463) submitter-processed matrix and platform annotation. GEO describes 30 idiopathic PAH and 41 healthy PBMC samples among 140 total. We use those 71 samples; the other 69 diagnoses are excluded from this contrast. The folder name `data/processed` is historical; these two GEO files are *inputs* to this pipeline. The paper PDF in `data/from paper/` is for local reference and is excluded from Git by `.gitignore`; none of the scripts reads it.

## Read the results

**Muscle, observed.** The nine PAH paper proteins were lower in a small resting muscle comparison. In separate MoTrPAC muscle RNA, positive endurance-versus-control estimates occur for 0/9 genes at 15–45 minutes, 7/9 at 3.5–4 hours (none BH q < 0.05), and 9/9 at 24 hours (seven BH q < 0.05). The 24-hour emphasis was chosen **after** inspecting the time course. Five corresponding genes belong to the GO oxidative-phosphorylation set. Its definition has 148 symbols; MoTrPAC's pathway result reports 139 measured genes. That RNA test is down early, inconclusive at 3.5–4 hours, and up at 24 hours. None of the nine matching MoTrPAC proteins has BH q < 0.05 at any of the three sampled protein times (0/27 matches).

**Muscle, interpretation and hypothesis.** The delayed healthy-cohort RNA pattern aligns with a disease-linked mitochondrial hypothesis. RNA abundance is not the rate of transcription, and higher RNA does not establish higher protein, improved muscle function, or the same response in patients. A patient study would need paired PAH and comparison muscle RNA, protein, and function after exercise, with a disease-by-exercise test.

**Blood, observed.** Script 07 fits `IPAH − Healthy` in PBMCs. Positive logFC/t means higher signal in IPAH. It tests 48,803 probes, maps 29,255 unambiguous probes, and selects one probe per 19,593 gene symbols. Script 08 matches 13,055 genes to MoTrPAC **whole-blood** RNA at each of six times. Spearman rho is small: −0.079 and −0.029 during exercise, then +0.060 to +0.125 afterward. Negative means opposing ranks; positive means same-direction ranks. The 10,000 gene-label shuffle P values break gene-gene dependence and are exploratory. Script 09 tests 5,183 GO sets with 10–500 PAH genes; 63 have BH FDR < 0.05. Script 11 joins 3,358 GO labels per time. Of 40 set-time rows with BH FDR < 0.05 in **both source analyses**, 27 point in the same direction and 13 in opposite directions.

**Blood, interpretation and hypothesis.** The PBMC and whole-blood cohorts differ in cell mixtures, platform, disease status, and sampling. GO terms share genes, so 40 rows do not represent 40 independent mechanisms. The probe with largest absolute moderated t was selected using the same disease outcome; its attached BH value is **probe-level BH, not formal gene-level FDR**. Script 12 checks whether the rho pattern changes when 24 matched genes from two blood-related GO sets are removed, or when the highest-mean-expression probe represents each gene. These checks do not remove general cell-composition differences. A future patient exercise study would need a direct disease-by-exercise comparison.

The concise claim for a team pitch is: **these are cross-cohort molecular relationships that motivate a test in PAH patients; they are not evidence of a PAH treatment effect.**
