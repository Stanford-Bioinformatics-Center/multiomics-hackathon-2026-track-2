# Metabolomics: PAH plasma metabolites in MoTrPAC context

This module starts from a public PAH metabolomics study, **ST000763** (systemic sclerosis with and without PAH, plasma at rest and at peak exercise), analyses it, and maps its metabolites onto the healthy MoTrPAC blood and muscle metabolomics. It compares separate studies. It does not test whether exercise treats PAH.

## Data prep

`data/workbench/ST000763/` holds the two files downloaded from the [Metabolomics Workbench](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST000763&StudyType=MS&ResultType=1) REST API on 25 September 2026: `factors.json` (sample labels) and `data.json` (named-metabolite values). The study page says blood was drawn "at rest and again at peak exercise for each human subject". There are **218 samples** and **673 metabolite features** (552 unique RefMet names) from four LC–MS runs. Values are "Peak area normalized" and are log2 transformed before any test.

Each sample carries four labels: `Exercise` (No = rest, Yes = peak), `Group`, `Sex` and `Type`. **`Type` matters.** Every PAH, borderline-pressure and normal-pressure sample is `Cath` (drawn in the catheterization setting). Every healthy and low-risk SSc sample is `Non-invasive`:

| Group | Setting | Rest | Peak | Inferred pairs |
|---|---|---:|---:|---:|
| Healthy | Non-invasive | 37 | 37 | 37 |
| Low-risk SSc | Non-invasive | 22 | 22 | 22 |
| Normal pressures (SSc) | Cath | 10 | 10 | 10 |
| Borderline pressures (SSc) | Cath | 21 | 21 | 21 |
| SSc-PAH | Cath | 20 | 18 | 18 |

So **PAH versus healthy always mixes three things: PAH, scleroderma and sampling setting.** The same-setting comparison is PAH versus normal-pressure SSc.

The Workbench files have **no participant IDs**. A rest sample is paired with the next-numbered peak sample from the same group, sex and setting. This gives **108 inferred pairs**; they are plausible, not verified.

`data/refmet/` holds the RefMet metabolite-class list copied from the `MotrpacHumanPreSuspensionAnalysis` 2.0.8 source. `data/source_manifest.json` records every input's SHA-256 hash; script 00 checks them.

MoTrPAC data come from the `MotrpacHumanPreSuspensionAnalysis` **2.0.8** R package (data collection c2.0). Three of its precomputed endurance contrasts are used:

| Contrast | Meaning | Use here |
|---|---|---|
| `EE-CON` | (exercise post − pre) − (control post − pre) | the exercise effect |
| `EE-EE` | exercise post − pre only | like-for-like with patient peak − rest |
| `CON-CON` | control post − pre only | drift without exercise (fasting, time of day, handling) |

## Scripts

**PAH data: what differs in PAH?**

0. `00_verify_inputs.py` checks the input hashes and that there are 218 samples and 673 features.
1. `01_prepare_st000763.py` reads the labels and values, builds the **108 inferred rest/peak pairs**, and writes the sample sheet (`01_st000763_samples.csv`), the pair list (`01_inferred_pairs.csv`) and the group-by-setting table above (`01_group_setting_counts.csv`). It saves the log2 matrix and feature table to `data/processed/` for the next scripts.
2. `02_compare_exercise_response.py` asks whether **exercise changes a metabolite differently in PAH**. For each pair, change = log2(peak) − log2(rest). It compares mean changes between groups, e.g. (PAH peak − rest) − (healthy peak − rest), with a Welch t test and BH correction. Result: **560 eligible metabolites, 0 with q < 0.05** (smallest q 0.84). Only 17 have raw P < 0.05, fewer than the ~28 expected by chance. PAH versus normal-pressure, borderline and all SSc without PH also give 0. Within PAH alone, four metabolites change at q < 0.05: **succinic acid rises** (+0.75 log2) and three long-chain fatty acids fall. Succinate also rises in healthy participants, so it is an exercise response, not a PAH one. No MoTrPAC data are used here.
3. `03_compare_resting_levels.py` asks whether **resting levels differ**, using rest samples only (no pairing needed):
   - **PAH vs healthy: 43 of 574 features at q < 0.05** (41 metabolites): 26 higher in PAH (mostly free fatty acids, about 2–3-fold) and 17 lower (triglycerides, phosphatidylcholines). The result survives rank-based testing (69), per-sample centering (45) and sex adjustment (37).
   - **PAH vs normal-pressure SSc (same Cath setting): 0.** Normal → borderline → PAH trend within Cath: **0.**
   - Normal-pressure SSc already carries a median **90%** of each PAH-versus-healthy difference (41 of 43 carry at least half).
   - Class level: saturated and unsaturated fatty acids are higher in PAH than healthy.

   *Question: is the PAH-versus-healthy difference PAH, scleroderma, or where the blood was drawn?* The data point to setting and/or scleroderma. A possible, **unverified** cause is heparin given during catheterization, which releases free fatty acids from triglycerides.

**MoTrPAC context: how do these metabolites behave in healthy people?**

4. `04_export_motrpac_metabolomics.R` exports MoTrPAC's precomputed `EE-CON`, `EE-EE` and `CON-CON` results for blood (**1,143 features × 6 times**) and muscle (**891 features × 3 times**): **28,593 rows**, each with `logFC` and BH-adjusted p value. No model is fitted. Output: `data/processed/motrpac_metabolite_contrasts.csv.gz` and `04_motrpac_export_counts.csv`.
5. `05_align_exercise_with_motrpac.py` asks how **similar each group's rest-to-peak change is to MoTrPAC's change** during exercise (20 and 40 minutes, closest to "peak"). Metabolites are matched by exact name: MoTrPAC's `refmet_name`, or its `feature_id` when the RefMet name is blank (mostly lipids named like "PC 40:4"). This matches **340 names** (178 through `feature_id`). Spearman rho across matched metabolites is descriptive only.

   | MoTrPAC contrast, 20 min | PAH rho (all 339) | Healthy rho (all 339) | PAH rho (RefMet-only 161) | Healthy rho (RefMet-only 161) |
   |---|---:|---:|---:|---:|
   | `EE-CON` | 0.24 | 0.12 | 0.44 | 0.27 |
   | `EE-EE` (like-for-like) | 0.08 | 0.19 | 0.29 | 0.30 |
   | `CON-CON` | −0.22 | −0.04 | −0.18 | −0.10 |

   With the like-for-like `EE-EE` contrast, **PAH is no more MoTrPAC-like than healthy**. PAH's apparent extra likeness under `EE-CON` comes from the control side: PAH changes run slightly opposite to MoTrPAC's resting-control drift, and subtracting that drift inflates the match. The RefMet-only row reproduces the earlier 162-name analysis (PAH 0.44 → 0.29).
6. `06_build_context_table.py` adds MoTrPAC context to **any PAH metabolite list** in one standard format (see "Adding a paper's metabolite list" below). It already includes script 03's resting PAH-versus-healthy list. For each metabolite it reports: the MoTrPAC exercise effect at each time, peak time, whether it is still changed at 3.5–24 h, the largest control drift, muscle response, and a label. A label needs **BH q < 0.05 and |log2 FC| ≥ 0.5** (about 1.4-fold); with q alone, 78% of all metabolites "drift" in MoTrPAC controls, which says nothing.
   - 41 resting hits: **27 match MoTrPAC blood; 24 are "stable in MoTrPAC blood"**, 3 triglycerides "drift without exercise", 14 have no MoTrPAC match (e.g. arachidonic acid is measured only in MoTrPAC muscle).
   - The hits are **no more exercise-sensitive than the study's other matched metabolites** (0% vs 4.1%, Fisher P 0.61) or drift-sensitive (11% vs 6.8%, P 0.42).
   - PAH-versus-healthy differences are a median **3.1 times** (range 1.1–16) MoTrPAC's largest control drift. Free fatty acids are mostly ≥ 2.5 times; the three triglycerides only 1.3–1.6 times, so fasting or time of day could contribute to those.

   Outputs: `06_context_table_hits.csv` (the readable table), `06_context_table_all.csv` (hits plus background) and `06_context_summary.json`.
7. `07_plot_resting_setting_effect.py` plots the 12 strongest resting results in every group, relative to healthy (`07_resting_setting_effect.png` / `.svg`, data in `07_resting_setting_effect_data.csv`). The three Cath groups sit together far from healthy; low-risk SSc sits near healthy. PAH does not stand apart from normal-pressure SSc.

## Adding a paper's metabolite list

Save a CSV in `data/pah_metabolite_lists/` with these columns, then rerun script 06:

| Column | Meaning |
|---|---|
| `list_id` | short ID for the list, e.g. `coursen2024_ssc_pah` |
| `source` | paper citation or dataset |
| `source_type` | `paper_reported` or `our_reanalysis` |
| `comparison` | e.g. SSc-PAH vs SSc without PH, rest |
| `metabolite_name` | name as the paper reports it |
| `refmet_name` | RefMet name (use the [RefMet name converter](https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_form.php)); blank uses `metabolite_name` |
| `pah_log2_effect` | PAH minus comparator on the log2 scale; blank if not reported |
| `pah_direction` | `up` or `down` |
| `pah_q_value` | blank if not reported |
| `is_reported_hit` | `True` for reported PAH-altered metabolites; `False` for other measured metabolites (background) |
| `caveat` | design limits, e.g. sampling setting |

Papers usually list only their significant metabolites. Without `False` background rows, the table still annotates each marker, but it cannot say whether PAH markers are more exercise-sensitive than metabolites in general.

## Run it

Use Python 3.10+ with the packages in `requirements.txt` (`python -m pip install -r requirements.txt`), and R 4.4+ with `MotrpacHumanPreSuspensionAnalysis` 2.0.8 for step 04. From any folder:

```sh
python "MoTrPAC Hackathon/Metabolomics/run_pipeline.py"
python "MoTrPAC Hackathon/Metabolomics/run_pipeline.py" --skip-r              # reuse the bundled MoTrPAC export
python "MoTrPAC Hackathon/Metabolomics/run_pipeline.py" --r-library PATH      # package in its own library
```

The local verification used Python 3.12.0, R 4.4.2 and `python run_pipeline.py --r-library ../../reference_code/R-library-c2`. Each numbered script also runs on its own after the steps before it. `--skip-r` lets you rerun everything without the MoTrPAC R package, because step 04's export is saved in `data/processed/`.

## Concise conclusion

In ST000763, no metabolite's acute exercise response differs between SSc-PAH and healthy people or other SSc groups. At rest, PAH differs strongly from healthy people, but the normal-pressure SSc group sampled in the same catheterization setting shows almost the same differences, and PAH versus that group gives none. MoTrPAC shows that most of these metabolites are stable in healthy people during exercise and without it, and that the resting differences are several times larger than normal physiological drift. So recent exercise or time of day does not explain them; sampling setting and/or scleroderma is the leading explanation, not PAH. The practical message for PAH metabolite biomarkers is to match sampling conditions and record recent activity. These are cross-study descriptions from inferred pairs and public labels, not proof of mechanism.
