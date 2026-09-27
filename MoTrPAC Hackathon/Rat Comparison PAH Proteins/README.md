# Rat Comparison PAH Proteins

This folder starts with the nine proteins in the **Downregulated proteins** part of Malenfant et al.'s Table 2. It maps each human paper entry to one measured MoTrPAC rat gastrocnemius protein, checks the rat sample groups, fits sex-adjusted training comparisons, and produces the nine-panel training-course plot. It uses **protein abundance only**. The plot describes healthy rats; it does not measure an exercise response in people with PAH.

## Data prep

`data/paper/pah_lower_proteins_malenfant2015.csv` transcribes the nine downregulated rows from Table 2 of Malenfant et al., including the paper's UniProt accession, symbol, PAH/control ratio, and reported P value. Table 2 is on page 5 of the local article PDF used to prepare this CSV. The publisher PDF is **not required or distributed with the shared pipeline**; use the article DOI under Sources if you need to check the table. For example, NDUFA9 has a ratio of **0.71**: measured abundance in the resting PAH group was about **29% lower** than in the paper's healthy controls. A ratio below 1 identifies a downregulated paper entry; it does not predict the rat exercise response.

`data/rat/` contains the official MoTrPAC `MotrpacRatTraining6moData` version 2.0.0 objects used here. `source_manifest.json` records the source repository, pinned commit, and SHA-256 hash for each `.rda` file. `PROT_SKMGN_NORM_DATA` is sample-level normalized gastrocnemius protein abundance. `PHENO` describes the rats. `RAT_TO_HUMAN_GENE` and `FEATURE_TO_GENE_FILT` provide ortholog and feature mappings. `PROT_SKMGN_DA` contains MoTrPAC's published sex-specific differential-abundance results. The data are copied here so this folder can run without depending on `module_staging`.

## Run the pipeline

Install **Python 3.10 or newer** and **R 4.4 or newer**. Only Python's standard library and base R are used. From a clone of the repository, change into this module folder and run the platform-independent runner:

```sh
cd "MoTrPAC Hackathon/Rat Comparison PAH Proteins"
python scripts/run_all.py
```

The runner checks the CSV and the SHA-256 hashes of all rat data files, finds `Rscript` on `PATH` (or in a standard Windows R installation), and runs Steps 01–05. If Rscript is elsewhere, use `python scripts/run_all.py --rscript "/path/to/Rscript"`. The included rat data allow this to run without downloading new data or installing an MoTrPAC R package.

To run **one step at a time**, stay in this module folder and run the commands below. `01` must finish before `02`, and so on. If `Rscript` is not on `PATH`, use its full executable path in place of `Rscript`. On Windows, `scripts/run_all.ps1` is also available; it accepts `-RscriptPath`.

```sh
python scripts/run_all.py --verify-only
Rscript scripts/01_map_paper_to_rat.R
Rscript scripts/02_prepare_rat_samples.R
Rscript scripts/03_estimate_training_course.R
Rscript scripts/04_compare_pooling_choices.R
Rscript scripts/05_plot_training_course.R
```

## Scripts and outputs

1. **`run_all.py --verify-only`** checks the nine paper rows and verifies every rat `.rda` file against the manifest hash. The Windows alternative is `00_verify_inputs.ps1`. Stop if a hash fails; the following results would no longer have the stated provenance.
2. **`01_map_paper_to_rat.R`** makes `output/01_paper_to_rat_mapping.csv`. It retains each paper UniProt accession for audit, applies an explicit paper-symbol-to-current-human-symbol crosswalk, then uses MoTrPAC's human-to-rat ortholog table and rat feature table. ATP5L is now **ATP5MG**; ATP5B is now **ATP5F1B**. The script demands one rat ortholog and one measured feature for each paper protein. It also exports `01_official_sex_specific_weekly.csv`: nine proteins × two sexes × four weeks = **72 published contrast rows**. `logFC` is trained minus sedentary. `adj_p_value` is MoTrPAC's published adjusted P value; `selection_fdr` tests a different training-selection question and should not be called an 8-week P value.
3. **`02_prepare_rat_samples.R`** matches each selected protein value to the rat's sex and training group. `02_rat_protein_samples.csv` contains nine proteins × 57 rats = **513 rows**; `02_sample_counts.csv` shows the group sizes. There are 47 trained protein samples across weeks 1, 2, 4, and 8, and 10 sedentary samples. The sedentary muscle samples were collected at the **8-week endpoint only**.
4. **`03_estimate_training_course.R`** fits `normalized log2 protein ~ sex + group` once per protein. It saves 45 sex-adjusted group means and 95% model confidence intervals in `03_sex_adjusted_group_means.csv`. Equal weighting of female and male predictions stops unequal group counts from changing what “sex-combined mean” means. `03_training_vs_sedentary_by_week.csv` has the trained-minus-sedentary estimates, P values, and exploratory BH q values across the nine chosen proteins **within each week**. `03_training_duration_patterns.csv` checks whether the four trained-group means strictly rise or fall and fits a separate, exploratory linear week slope among trained rats. A significant slope is not the same as every successive mean increasing.
5. **`04_compare_pooling_choices.R`** answers the earlier trained-versus-untrained question in `04_training_comparisons.csv`. It reports all trained weeks combined versus sedentary, then repeats at the 8-week endpoint with females alone, males alone, and both sexes with sex as a model covariate. Each comparison adjusts its nine target P values with BH. The 8-week version holds the sampling endpoint constant. The all-week version averages training durations and compares some earlier trained rats with later sedentary rats.
6. **`05_plot_training_course.R`** reads the Step 03 means and writes `05_rat_protein_training_course.png`. The magenta line joins the trained groups at weeks 1, 2, 4, and 8; each point comes from **different rats**. The blue dashed line and pale band repeat the **one 8-week sedentary estimate and its 95% confidence interval** as a visual reference. The blue dot at week 8 is the only sampled sedentary time. There are no sedentary measurements at weeks 1, 2, or 4.

## Reading the results

At 8 weeks, all nine sex-adjusted trained-minus-sedentary estimates are positive; seven have **targeted BH q < 0.05**. In the sex-combined trained means, **five of nine** rise at every successive sampled duration. This is a pattern across different animals, not a within-rat progression. The trained-only linear trend model is exploratory and assumes a straight-line duration effect.

Pooling sexes at the same 8-week endpoint detects **7/9**, versus **6/9 female-only** and **4/9 male-only**, using the same targeted nine-test BH method. Pooling may improve precision when responses align, but does not increase the underlying effect size. In contrast, the official MoTrPAC sex-specific, much broader BY-adjusted analysis detects **3/9 female** and **0/9 male** at week 8. These two adjustment schemes answer different multiple-testing questions and their hit counts should not be compared as if pooling were the only change.

The all-week pooled comparison detects **4/9** with targeted BH q < 0.05. It mixes training durations, so use the time-course plot and the 8-week comparison when discussing longer training. Neither analysis establishes that exercise reverses protein changes in PAH patients or treats PAH. Malenfant's PAH samples and MoTrPAC's healthy-rat training samples come from different species, cohorts, and study designs.

## Sources

- Malenfant et al., *Journal of Molecular Medicine* (2015), Table 2, DOI [10.1007/s00109-014-1244-0](https://doi.org/10.1007/s00109-014-1244-0).
- MoTrPAC Study Group, *Nature* 629, 174–183 (2024), DOI [10.1038/s41586-023-06877-w](https://doi.org/10.1038/s41586-023-06877-w); rat data package [MotrpacRatTraining6moData](https://github.com/MoTrPAC/MotrpacRatTraining6moData), commit recorded in `data/rat/source_manifest.json`.
