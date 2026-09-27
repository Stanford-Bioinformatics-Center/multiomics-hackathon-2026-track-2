# Robustness matrix report

`tests/matrix/run_matrix.py`: 52 cases, run with `mprobe run` (1,000 random sets per null) in parallel on the desktop node. Assertions: `tests/test_matrix.py`. Metabolite cases: `tests/test_matrix_metab.py`.

## Summary

- Exit 0: 52 / 52.
- Reports that validate (no broken refs, captions present, no empty tables): 52 / 52.
- Runtime per run: median 5.8 s, max 11.0 s (limit 90 s).

## By group

| group | cases | median_runtime_s | median_class_pct | cases_with_setlevel_verdicts |
|---|---|---|---|---|
| edge cases | 2 | 4.8 | 71.7 | 0 |
| example: disease | 6 | 7.2 | 71.8 | 3 |
| example: negative control | 1 | 4.8 | 61.7 | 0 |
| example: pathway-class control | 1 | 5.1 | 66.1 | 1 |
| example: positive control | 1 | 6.7 | 2.5 | 1 |
| id types | 9 | 5.5 | 83.3 | 0 |
| seed | 2 | 5.4 | 64.0 | 1 |
| size x composition | 28 | 6.5 | 50.7 | 2 |
| weights | 2 | 5.8 | 83.3 | 0 |

## All cases

n_opposed_fdr / n_same_fdr = core comparisons (of 52) with cameraPR BH FDR < 0.05 opposed / same direction; median_pct_* = median null percentile over the 52; rat_gn_prot_f8_t = signed cameraPR t in rat gastrocnemius protein, female 8 wk; n_flips = sensitivity-table conclusions that change under an alternative setting.

| case | exit | runtime_s | n_rows | n_mapped | n_unmapped | n_columns_no_genes | n_opposed_fdr | n_same_fdr | median_pct_abund | median_pct_class | rat_gn_prot_f8_t | n_flips |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| size5_random_half | 0 | 5.30 | 5 | 5 | 0 | 0 | 0 | 0 | 56.70 | 54.20 | 0.99 | 4 |
| size5_random_up | 0 | 5.40 | 5 | 5 | 0 | 0 | 0 | 0 | 89.35 | 86.30 | 0.68 | 1 |
| size5_random_down | 0 | 4.30 | 5 | 5 | 0 | 0 | 0 | 0 | 10.65 | 13.70 | -0.95 | 2 |
| size25_random_half | 0 | 6.30 | 25 | 25 | 0 | 0 | 0 | 0 | 57.05 | 58.90 | -0.14 | 0 |
| size25_random_up | 0 | 6.30 | 25 | 25 | 0 | 0 | 0 | 1 | 40.85 | 54.00 | -2.12 | 20 |
| size25_random_down | 0 | 6.50 | 25 | 25 | 0 | 0 | 1 | 0 | 59.15 | 46.00 | 1.72 | 3 |
| size100_random_half | 0 | 8.20 | 100 | 100 | 0 | 0 | 0 | 0 | 56.30 | 50.20 | 0.44 | 1 |
| size100_random_up | 0 | 8.00 | 100 | 100 | 0 | 0 | 0 | 0 | 53.20 | 41.65 | -0.26 | 10 |
| size100_random_down | 0 | 8.00 | 100 | 100 | 0 | 0 | 0 | 0 | 46.80 | 58.35 | -0.46 | 0 |
| size500_random_half | 0 | 10.20 | 500 | 500 | 0 | 0 | 0 | 0 | 47.30 | 42.50 | -0.54 | 0 |
| size500_random_up | 0 | 8.40 | 500 | 500 | 0 | 0 | 0 | 0 | 57.20 | 58.65 | -1.23 | 39 |
| size500_random_down | 0 | 9.90 | 500 | 500 | 0 | 0 | 0 | 0 | 42.80 | 41.35 | 0.15 | 1 |
| size5_mito_half | 0 | 4.90 | 5 | 5 | 0 | 0 | 0 | 0 | 34.25 | 45.40 | -1.25 | 0 |
| size25_mito_half | 0 | 5.50 | 25 | 25 | 0 | 0 | 0 | 0 | 34.20 | 38.25 | -0.40 | 13 |
| size100_mito_half | 0 | 6.60 | 100 | 100 | 0 | 0 | 0 | 0 | 44.10 | 48.00 | -1.36 | 32 |
| size500_mito_half | 0 | 9.90 | 500 | 500 | 0 | 0 | 0 | 0 | 74.75 | 53.90 | 0.48 | 26 |
| size5_secreted_half | 0 | 4.30 | 5 | 5 | 0 | 0 | 0 | 0 | 48.50 | 51.25 | -0.96 | 0 |
| size25_secreted_half | 0 | 5.70 | 25 | 25 | 0 | 0 | 0 | 0 | 55.30 | 54.55 | -0.83 | 0 |
| size100_secreted_half | 0 | 6.50 | 100 | 100 | 0 | 0 | 0 | 0 | 74.90 | 77.15 | 0.49 | 3 |
| size500_secreted_half | 0 | 11.00 | 489 | 489 | 0 | 0 | 0 | 0 | 43.60 | 56.15 | -0.58 | 11 |
| size5_contractile_half | 0 | 4.40 | 5 | 5 | 0 | 0 | 0 | 0 | 34.25 | 33.05 | -1.61 | 0 |
| size25_contractile_half | 0 | 5.40 | 25 | 25 | 0 | 0 | 0 | 0 | 61.35 | 61.00 | -0.01 | 6 |
| size100_contractile_half | 0 | 6.70 | 100 | 100 | 0 | 0 | 0 | 0 | 62.65 | 59.65 | -0.35 | 32 |
| size500_contractile_half | 0 | 7.50 | 236 | 236 | 0 | 0 | 0 | 0 | 40.90 | 39.30 | -0.64 | 41 |
| size5_ribosomal_half | 0 | 4.90 | 5 | 5 | 0 | 0 | 0 | 0 | 44.15 | 50.05 | -0.25 | 3 |
| size25_ribosomal_half | 0 | 5.40 | 25 | 25 | 0 | 0 | 0 | 0 | 55.45 | 45.10 | -0.11 | 5 |
| size100_ribosomal_half | 0 | 6.60 | 100 | 100 | 0 | 0 | 0 | 0 | 68.35 | 58.65 | -0.04 | 3 |
| size500_ribosomal_half | 0 | 7.80 | 220 | 220 | 0 | 0 | 0 | 0 | 36.45 | 38.85 | -0.99 | 28 |
| id_human_symbols | 0 | 5.30 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| id_uniprot_column | 0 | 5.20 | 25 | 16 | 9 | 0 | 0 | 0 | 74.25 | 82.90 | 0.87 | 1 |
| id_uniprot_in_symbol | 0 | 5.60 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| id_ensembl | 0 | 5.50 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| id_rat_symbols | 0 | 5.30 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| id_mixed_case | 0 | 5.70 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| id_duplicates | 0 | 5.80 | 33 | 22 | 11 | 0 | 0 | 0 | 80.05 | 82.55 | 1.61 | 1 |
| id_20pct_unknown | 0 | 5.20 | 30 | 24 | 6 | 0 | 0 | 0 | 80.40 | 84.55 | 0.87 | 1 |
| absent_from_proteomics | 0 | 5.60 | 25 | 25 | 0 | 6 | 0 | 0 | 74.10 | 71.65 | 1.87 | 1 |
| empty_overlap_tissue | 0 | 4.00 | 8 | 8 | 0 | 46 | 0 | 0 |  |  |  | 0 |
| weights_present | 0 | 5.80 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| weights_absent | 0 | 5.90 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| direction_words | 0 | 5.50 | 25 | 25 | 0 | 0 | 0 | 0 | 79.45 | 83.30 | 0.87 | 2 |
| example_pah_muscle_malenfant2015 | 0 | 7.30 | 25 | 25 | 0 | 0 | 19 | 2 | 90.45 | 64.55 | 4.62 | 24 |
| example_pah_blood_cheadle2012_eds | 0 | 4.90 | 22 | 22 | 0 | 0 | 0 | 1 | 37.00 | 48.25 | 0.90 | 25 |
| example_hostrup2022_hiit_proteome | 0 | 6.70 | 125 | 123 | 2 | 0 | 2 | 36 | 0.00 | 2.45 | -9.11 | 11 |
| example_random_matched | 0 | 4.80 | 25 | 25 | 0 | 0 | 0 | 0 | 68.75 | 61.70 | 0.97 | 3 |
| example_random_mito9 | 0 | 5.10 | 9 | 9 | 0 | 0 | 4 | 0 | 94.05 | 66.05 | 2.76 | 11 |
| example_type2_diabetes_muscle_mootha2003 | 0 | 5.50 | 87 | 87 | 0 | 0 | 33 | 8 | 99.90 | 92.20 | 11.62 | 17 |
| example_heart_failure_lv_hannenhalli2006 | 0 | 8.60 | 600 | 571 | 29 | 0 | 0 | 0 | 98.15 | 98.80 | 0.53 | 28 |
| example_aging_muscle_liu2013 | 0 | 8.00 | 600 | 522 | 78 | 0 | 0 | 0 | 89.60 | 79.00 | 0.19 | 2 |
| example_septic_shock_blood_cvijanovich2008 | 0 | 7.20 | 600 | 510 | 90 | 0 | 0 | 0 | 43.75 | 30.45 | -0.85 | 31 |
| seed2_pah_muscle_malenfant2015 | 0 | 6.10 | 25 | 25 | 0 | 0 | 19 | 2 | 89.35 | 65.30 | 4.62 | 23 |
| seed2_random_matched | 0 | 4.70 | 25 | 25 | 0 | 0 | 0 | 0 | 65.80 | 62.65 | 0.97 | 4 |

## Seed reproducibility (seed 20260926 vs 12345)

Non-null quantities must be identical; null percentiles vary with the random draws.

| signature | quantity | kind | max_abs_diff | median_abs_diff |
|---|---|---|---|---|
| pah_muscle_malenfant2015 | camera_t | non-null | 0.00 | 0.00 |
| pah_muscle_malenfant2015 | n_opposed | non-null | 0.00 | 0.00 |
| pah_muscle_malenfant2015 | sign_p | non-null | 0.00 | 0.00 |
| pah_muscle_malenfant2015 | pct_abund | null-based | 5.10 | 0.50 |
| pah_muscle_malenfant2015 | pct_class | null-based | 4.80 | 1.05 |
| random_matched | camera_t | non-null | 0.00 | 0.00 |
| random_matched | n_opposed | non-null | 0.00 | 0.00 |
| random_matched | sign_p | non-null | 0.00 | 0.00 |
| random_matched | pct_abund | null-based | 6.80 | 1.05 |
| random_matched | pct_class | null-based | 5.80 | 1.35 |

## Failures

None.

## Found and fixed

- **Crash when RNA–protein ρ is exactly ±1** (cases `size5_secreted_half`, `absent_from_proteomics`, first matrix
  run): the Fisher-z CI bound is clipped at ±0.999999, so the error bar went slightly negative and matplotlib raised.
  Fixed in `figures.discordance` and `discord.fig_agreement` (error bars clipped at 0).
- **Empty tables when no ranking unit has enough measured genes** (`empty_overlap_tissue`): section F rendered
  empty tables, which the validator rejects. Fixed in `run.write`: an empty table is replaced by a one-row
  "no rows: …" table. The headline shows "no measured genes" cells, not an error.
- **Negative control above the class null** (found by the Phase 1 suite, confirmed here): the class null matched
  only on MitoCarta / complex / secreted flags, so a high-abundance random set sat above the 98th percentile in
  6 of 52 columns. The class null now also matches on abundance tertile (`nulls.column_null`).

