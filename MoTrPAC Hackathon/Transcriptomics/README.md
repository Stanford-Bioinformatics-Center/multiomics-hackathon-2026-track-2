# BioCarta pathways in MoTrPAC blood RNA

This module follows the six BioCarta pathways listed as lower in patient blood cells in Cheadle et al.'s Table 2, finds the same **pathway names** in MoTrPAC's healthy-volunteer blood RNA results, and makes the bubble plot in `outputs/biocarta_all_six_pathways.png`. It compares two separate studies. It does not test whether exercise changes these pathways in PAH patients.

## Data prep

`data/paper/Cheadle2012.pdf` is the paper (*PLoS ONE* 7:e34951, 2012; [DOI 10.1371/journal.pone.0034951](https://doi.org/10.1371/journal.pone.0034951)). `data/paper/table2_selected_pathways.csv` is a checked **manual transcription** of its eight selected Table 2 rows. The paper calls them “selected pathway results” and says a negative PAGE score means lower average expression in the disease group than in healthy controls. For example, the TOB1 pathway's IPAH versus control PAGE Z-score is −4.620. That number describes the paper's patient comparison; it is not an MoTrPAC exercise effect.

Six Table 2 rows have the annotation `BioCarta`. The other two use `GEArray` and `SABiosciences`. Selecting BioCarta gives six exact names in MoTrPAC's `C2/BIOCARTA` collection. This does **not** establish that the paper and MoTrPAC used identical gene members, that the other two have no biological overlap with MoTrPAC, or why the colleague chose BioCarta. The supplied example plot showed five rows, but the reason it omitted ARE/NRF2 is unknown. This module includes all six, as requested.

`data/motrpac/CAMERA_RESULTS.rda` and `data/motrpac/DESCRIPTION` are copied from the project's official `MotrpacHumanPreSuspensionAnalysis` version **2.0.8** source. The `.rda` contains already computed pathway tests. The source file hashes are in `data/source_manifest.json`; script 00 checks them before proceeding. No download or patient-level MoTrPAC reanalysis is needed.

## Scripts

0. `scripts/00_select_paper_biocarta.py` checks the four source-file hashes and the eight manually transcribed rows. It selects the six rows explicitly labelled BioCarta and writes `outputs/paper_table2_biocarta.csv`. This makes the selection visible instead of hiding it in the figure code.
1. `scripts/01_export_motrpac_biocarta.R` reads MoTrPAC's **precomputed** CAMERA results. It keeps blood `transcript-rna-seq`, exercise versus time-matched control, the `C2/BIOCARTA` collection, and the six selected names. It writes `outputs/motrpac_biocarta_camera.csv`: **60 rows = six sets × (six endurance exercise times + four resistance exercise times)**. The script does not run CAMERA again.
2. `scripts/02_map_paper_to_motrpac.py` joins the paper list to MoTrPAC by the exact BioCarta set name. It checks that all 60 set-by-time cells are present and writes `outputs/paper_motrpac_biocarta_mapped.csv` plus `outputs/pathway_mapping_qc.csv`. The four paper PAGE scores travel with each corresponding MoTrPAC result for inspection; the script does not subtract or correlate their magnitudes.
3. `scripts/03_plot_biocarta_bubbles.py` makes `outputs/biocarta_all_six_pathways.png` and `.svg`. Its rows are clustered by their ten MoTrPAC CAMERA Z-scores. Its columns remain in time order, endurance first and resistance second. Dot color is signed MoTrPAC CAMERA Z; dot **area** increases as the original MoTrPAC BH-adjusted P value decreases. This color/size choice reconstructs the supplied image provisionally because its full caption and plotting code were unavailable.

## Run it

In PowerShell, from any folder:

```powershell
& 'C:\Users\vanes\Downloads\MoTrPAC\MoTrPAC Hackathon\Transcriptomics\run_pipeline.ps1'
```

The runner uses Rscript (tested with R 4.4.2) and `python` (tested with Python 3.12.0). It needs the Python packages in `requirements.txt`. If one is missing, install them in your chosen Python environment with `python -m pip install -r requirements.txt` from this folder, then rerun. The R step uses base R only. Every script finds the module by its own location, so you can run each numbered file separately after its preceding step.

## Reading the figure

The paper's negative PAGE values mean lower pathway expression in the disease group than in controls across four **resting patient** comparisons. The MoTrPAC CAMERA Z values compare acute exercise with a matched nonexercise control in **healthy volunteers**. Red is an upward pathway tendency after exercise; blue is downward. A larger dot means a smaller MoTrPAC BH-adjusted P value; the q values were adjusted in MoTrPAC's original testing family, not recalculated for these six selected sets.

For example, the paper reports TOB1 at PAGE Z = −4.620 in IPAH versus controls. Its MoTrPAC dots include both positive and negative CAMERA Z values across time. The sign comparison is descriptive because PAGE and CAMERA use different methods, samples, and gene-set versions. A positive pathway Z also does not mean every member gene increased.

**Concise conclusion:** Among the 60 selected MoTrPAC pathway–time cells, 49 CAMERA Z-scores are positive and 11 negative; 21 positive cells have original BH q < 0.05, and no negative cell does. The paper reports lower pathway expression in patient PBMCs. This cross-study pattern generates a question about whether patients have different acute responses. It is not a disease-by-exercise test, evidence of immune restoration, or evidence that exercise treats PAH. PBMCs versus whole blood, blood-cell composition, shared genes among pathways, and possibly different gene-set membership limit interpretation.
