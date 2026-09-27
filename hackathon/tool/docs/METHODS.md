# Methods and module details

The full method notes that used to open each module. Each source file now keeps a short summary and points here.

## compare.py

`mprobe compare`: several disease signatures side by side (the generality view).

Rows = signatures, columns = layer x time: human vastus lateralis 24 h RNA and protein (EE, RE), rat gastrocnemius
protein 4 and 8 wk and RNA 8 wk (F, M), and rat heart protein 8 wk (F, M) so heart signatures have their tissue.
Cell = signed cameraPR t (+ = exercise opposes the signature) with the abundance-matched null percentile. The table
adds the pathway-class null percentile and each signature's MitoCarta share. The only interpretive sentence is
selected by rules on these numbers (see `lesson`).

## config.py

Config files that mirror the CLI flags (see tool/config.example.yaml).

YAML is read with PyYAML when installed; otherwise a flat subset is parsed here (top-level `run:` / `discord:`
blocks of `key: value`, lists as `[a, b]` or `a, b`). JSON files work everywhere. Keys use the CLI flag names with
underscores (context_groups, pool_sets, ...). CLI flags given explicitly win over the file.

## core.py

Scoring primitives shared by `run`, `library` and `discord`.

Definitions follow NARRATIVE.md / scripts/05_pah_grid.py exactly:
  agreement  = sign(exercise logFC) x signature direction; -1 = exercise moves the gene AGAINST the signature
  cell value = agreement x min(|stat|, cap), cap 4
  sign test  = scripts/05_pah_grid.py::sign_test (two-sided binomial, H0 0.5; zero / NA agreements dropped)
  cameraPR   = scripts/05_pah_grid.py::camera_pr (Python port of limma::cameraPR, inter.gene.cor 0.01) on the full
               column ranking; for a signed signature the member statistics are multiplied by -direction, so a
               POSITIVE t means the set is OPPOSED by exercise (up-half and down-half are also tested alone).

## discord.py

`mprobe discord`: omic-discordance report for one species x tissue (generalises scripts 04/06/08/11 + lag/).

Panels (order fixed): a timescale, b agreement by timepoint, c agreement by gene property, d which layer
responded and the detection-power decomposition, e protein-without-RNA calls, f does early RNA predict later
protein, g sex without thresholds, h caveats. Every number shown is also written to tables/*.csv.

## discord_extra.py

Two extra panels for `mprobe discord` (rat training, tissues with RNA and PROT): pure functions that return
tables + matplotlib figures; the caller renders them.

protein_without_rna(S, tissue, ...)
    Summary-statistics re-implementation of hackathon/protonly (PLAN.md, scripts/11_protein_only.py): replicated
    protein training responders whose RNA response is bounded below the protein effect ("protein-only"), with a
    label-shuffle analogue, a concordance null, a parametric bootstrap and the mirror (RNA-only).
lag_panel(S, tissue, ...)
    Summary-statistics port of hackathon/scripts/07_lag.py: does RNA at week w add held-out information about the
    8-week protein response beyond protein, and beyond same-time (8-week) RNA?

Both return dict(tables, figures, numbers, captions, notes):
    tables   : dict[str, pandas.DataFrame]
    figures  : dict[str, matplotlib.figure.Figure]   (the caller saves / closes them)
    numbers  : dict[str, float | int]                 (headline scalars)
    captions : dict[str, (title, how_to_read)]        (one entry per table and figure key)
    notes    : list[str]                              (method statements, deviations, comparison to the originals)

Every value is computed from the store (per gene x column logFC, stat, baseline_expr). SE = |logFC / stat|
(DESeq2 Wald z for rat RNA, limma t for rat PROT). Reference numbers quoted in notes are copied from
protonly/RESULTS.md and lag/RESULTS.md (join_table_v2, raw logFC) and are labelled as such.

## everywhere.py

Section F: the signature in every MoTrPAC comparison column (scripts/10_everywhere.py, generalised).

Ranking unit = tissue x layer x time (rat: the two sexes of one training week pooled; human: one contrast), counts
over the counted signature genes; sign_p = two-sided binomial of opposed vs 0.5; sign_q_bh = BH over all units;
units with fewer than min_n measured cells are kept in the CSV but left out of the ranking (10_everywhere: 10).
Non-exercise reference comparisons (control time course, baseline group differences, EE-RE) are included as a
negative control. The per-column cameraPR t (signed; + = opposed) is added as the primary set-level score.

## layers.py

Coverage + detection-power flags (section A) and RNA-vs-protein layer discordance (section E).

Detection-power flags follow scripts/08_discordance_model.py and model/discordance_report.pdf, which found that
"which layer responded" is largely which assay could see the gene:
  RNA baseline      human: limma AveExpr (log2 CPM); rat: log2(sedentary mean normalized count + 1)
                    (reference_average_intensity of the DESeq2 table); median over the tissue's columns, tertile
                    among all genes measured in that tissue's RNA columns.
  PROT missingness  rat: numNAs of the protein feature (median over columns); human: not in the public DA table.
  GTEx v8 TPM       median TPM in the matched human tissue (CFDE dataset), tertile among the tissue's RNA genes.
  n_collapsed       features collapsed into the gene (max over columns); > 1 = isoforms / multi-mapping.
Per-gene layer class (scripts/08, fdr_bh < cutoff in each layer): concordant (both significant, same sign),
opposite (both, opposite sign), RNA_only, PROT_only, ns.

## legacy.py

Bridge to the analysis scripts in hackathon/scripts/.

Those scripts run their whole pipeline at module level, so they cannot be imported directly. This module
parses each script, keeps only the named top-level functions and constant assignments, and executes those
nodes in a fresh namespace. The code that runs is the scripts' own source (no copies), so a change in the
scripts propagates to the tool; `SOURCES` records which names come from where, for provenance.json.

## library.py

MoTrPAC as a signature library: every store column exported as UP / DOWN gene sets (GMT + JSON index), and a
query that ranks all columns against a user signature.

Set definition (build):
  UP   = genes with fdr_bh < cutoff and logFC > 0 in the column, ranked by |stat| descending, first `cap` kept;
  DOWN = genes with fdr_bh < cutoff and logFC < 0, same ranking and cap.
  All 422 store columns (every kind and layer) are exported. Sets with 0 genes are not written to the GMT but are
  listed in the JSON index with n_genes = 0.

Query (per column):
  The cameraPR / sign-test scores come from everywhere.run (signed cameraPR t, + = the column moves the signature
  genes AGAINST the signature). Set overlap uses the counted signature genes split by direction (sig_up, sig_down)
  and the column's library sets (UP, DOWN, after the FDR cutoff and the cap):
    n_overlap_same    = |sig_up & UP|   + |sig_down & DOWN|
    n_overlap_opposed = |sig_up & DOWN| + |sig_down & UP|
    denominator       = |(sig_up | sig_down) | (UP | DOWN)|
    jaccard_same      = n_overlap_same / denominator
    jaccard_opposed   = n_overlap_opposed / denominator
  (0 when the denominator is 0).

## metab.py

Metabolomics (METAB) as a third layer: store view, metabolite signatures and a 3-layer pathway panel.

Rows come from store/metab.parquet (store.build_metab; schema in store/SCHEMA.md "## METAB layer"). They have the
contrasts.parquet schema, but for layer == "METAB" `gene_symbol_human` holds the RefMet NAME of a metabolite.

Pathway panel: for each pathway in store/pathway_map.csv and each comparison column of one tissue, the cameraPR t
(legacy.camera_pr, inter-gene correlation 0.01) of the pathway's gene set (RNA, PROT columns) or metabolite set
(METAB columns) against every other feature measured in that column. t > 0 = the set moves up with exercise
relative to the rest of the column. The map is many-to-many: a pathway has a gene side and a metabolite side, each
tested within its own layer. Genes and metabolites are never matched one-to-one, and metabolite pools are not flux.

## nulls.py

Calibration nulls for the opposition score (section D).

(i)  abundance-matched: each counted signature gene measured in the column is replaced by a random non-signature
     gene from the same abundance decile of that column (NARRATIVE R1b: decile of the column's own baseline, i.e.
     human AveExpr / rat sedentary mean count; columns without a baseline use GTEx v8 median TPM of the matched
     human tissue; genes without either form their own stratum). Drawn with replacement, gene by gene, in the
     same order as scripts/05_pah_grid.py R1b, and given the signature's direction vector.
(ii) pathway-class-matched: same, but the stratum is the gene's class = (MitoCarta3.0 member, GO:CC complex
     subunit, secreted/extracellular GO:CC proxy; the three flags of scripts/08_discordance_model.py) crossed with
     the abundance tertile of (i), so it is also abundance-matched. Answers "is it your signature, or its pathway
     class?".
(iii) pool null (optional, NARRATIVE R1a): random sets of the same size drawn without replacement from a fixed
     gene pool (e.g. KEGG OXPHOS + TCA genes measured in both layers), one draw list per dataset.
Scores: signed cameraPR t (primary) and n opposed. Percentile = 100 x P(null <= observed) (as in R1a); empirical
p = (1 + #null >= observed) / (B + 1) (as in R1b).

## ranked.py

Full disease rankings (gene + t over the whole genome), generalised from the team's PAH blood pipeline
(`MoTrPAC Hackathon/scripts` 07.5–12 on main: GSE33463 IPAH-minus-healthy PBMC ranking vs MoTrPAC blood RNA).

Three analyses, any disease, any MoTrPAC column:
1. Rank association (their script 08): Spearman rho between the disease statistic and the exercise statistic over
   shared genes; two-sided gene-label permutation p (10,000 shuffles; gene-gene dependence is not preserved, so p
   is exploratory, as they state). Added: Fisher-z 95% CI, BH across columns, and a calibration against the
   non-exercise reference contrasts of the same tissue x layer (control time course, baseline and EE-RE
   differences): how unusual is rho for this ranking when no exercise contrast is involved?
   Sign convention: rho > 0 = exercise moves genes the SAME way as the disease; rho < 0 = opposed.
2. Pathway concordance (their scripts 09–11): GO:BP sets from MotrpacHumanPreSuspensionAnalysis
   MOLECULAR_SIGNATURES (sets with 10–500 genes in the universe), cameraPR on the disease ranking and on each
   exercise column; BH within each side; a set is concordant (same direction) or discordant when BH < 0.05 on
   both sides. Added: the exercise side is computed with the same validated cameraPR on the store (so any
   tissue / species works), and the Spearman correlation of set t values is reported next to the counts.
   `exercise_side="precomputed"` uses the package's CAMERA_RESULTS instead (human only), which is exactly the
   team pipeline and reproduces its 40 jointly significant set x time rows (27 same, 13 opposite) for PAH blood.
3. A directional signature derived from the ranking (BH < 0.05, at most 250 genes per direction by |t|) feeds
   every other section of the report.

## render.py

HTML report rendering for motrpac_probe.

Produces one self-contained HTML file per run (inline CSS, inline SVG
figures, HTML tables, optional inlined Vega-Lite for interactive items) and
a static gallery index. Colours and fonts for the HTML live in
``templates/style.css``; matplotlib colours live in the constants below.

Report context contract (``render_report(context, out_html)``)
--------------------------------------------------------------
Plain-text fields are HTML-escaped by the template. To pass inline HTML in a
plain-text field wrap it in ``markupsafe.Markup``. Fields documented as
"HTML" are inserted verbatim (they come from ``df_to_html`` / ``fig_to_svg``).

``title`` : str
    Page title (h1).
``subtitle`` : str, optional
    One line under the title.
``meta`` : list of (label, value)
    Compact key-value line under the title, e.g. signature name, n genes,
    run date, store hash, git SHA, command.
``headline`` : dict
    The first block after the title. Keys ``title`` (str), ``how_to_read``
    (str), ``table_html`` (HTML), ``sentences`` (list of str).
``first_screen`` : list of str, optional
    Short lines (for example one caveat) shown right under the headline.
``sections`` : list of dict
    Each ``{"id", "title", "intro": [str], "items": [...]}``. The table of
    contents is built from these. An item is one of::

        {"kind": "figure", "title", "how_to_read", "svg" (HTML), "png": "figures/x.png", "id"}
        {"kind": "table", "title", "caption", "html" (HTML), "csv": "tables/x.csv"}
        {"kind": "text", "html" (HTML)}
        {"kind": "interactive", "title", "how_to_read", "spec_json" (Vega-Lite JSON str),
         "id", "fallback_svg" (HTML)}
        {"kind": "details", "summary", "items": [...]}   # collapsible block

    Figures and interactive items without a title or how-to-read line render
    a red "MISSING CAPTION" marker, which ``validate_report`` reports.
``glossary`` : list of (term, definition), optional
``provenance_json`` : str, optional
    Pretty-printed JSON, shown in a collapsed block at the end.
``toggles`` : dict of section id -> bool, optional
    Sections mapped to False are not rendered (``"glossary"`` and
    ``"provenance"`` may also be toggled off).
``tool_version`` : str, optional
    Defaults to ``motrpac_probe.__version__``; shown in the footer.

Paths in ``png`` / ``csv`` are relative to the HTML file; the caller writes
those files. The vendored Vega JS (``vendor/``) is inlined only when at
least one interactive item is present.

## signature.py

Read a signature CSV and map it to the human gene symbols used in the store.

Input columns (case-insensitive): gene_symbol and/or uniprot / ensembl / rat_symbol, direction (+1 up in disease,
-1 down; 'up'/'down' accepted), optional group, weight, source. Mapping order per row, first hit wins:
  1. gene_symbol exactly as given, then upper-cased, against the store's human symbols;
  2. a gene_symbol value that looks like an Ensembl (ENSG/ENSRNOG) or UniProt accession is mapped as such;
  3. rat symbol -> human ortholog (RGD table used by scripts/03_join.py), case-insensitive;
  4. the uniprot column (MoTrPAC human proteomics feature map; isoform suffix dropped if needed);
  5. the ensembl column; 6. the rat_symbol column.
Duplicated genes keep the first row if directions agree and are dropped if they conflict. Every input row is
reported with its mapping status.

## store.py

Build and load the MoTrPAC contrast store (tool/store/*.parquet). Schema: tool/store/SCHEMA.md.

Rows = one gene per comparison column (dataset x tissue x layer x contrast), exactly as scripts/03_join.py and
scripts/10_everywhere.py define them:
  * columns present in data/join_table_v2.csv are copied from it verbatim (source = "join_table_v2");
  * every other column (19 rat tissues, human muscle/adipose/blood, all contrasts incl. non-exercise references,
    PHOSPHO) is collapsed from the raw DA tables with the `bh` and `collapse` functions lifted from
    10_everywhere.py (source = "alltissue_v1"): BH over all tested features before collapsing, rat RNA on/off
    filter, max-|stat| feature per human gene.

