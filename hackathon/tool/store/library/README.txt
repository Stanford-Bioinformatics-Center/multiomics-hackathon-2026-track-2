MoTrPAC signature library (motrpac_probe)
==========================================

Files
  motrpac_library.gmt          gene sets, standard GMT format (one set per line:
                                 name<TAB>description<TAB>gene1<TAB>gene2...). Genes are human
                                 HGNC symbols (rat genes mapped to human orthologs as in the store).
  motrpac_library_index.json   one JSON object per set, including empty sets (n_genes = 0, not in
                                 the GMT): set name, column_id, species, dataset, tissue, layer,
                                 contrast, label, kind, sex, time, direction, n_genes,
                                 n_significant_before_cap, cutoff, cap, provenance.

Definition
  Every MoTrPAC comparison column in the store (422 columns: rat training vs sedentary, human acute
  exercise vs control and post-pre, and the non-exercise reference comparisons; RNA, protein,
  phosphosite and Olink layers) gives two sets:
    UP   = genes with BH FDR < 0.05 and logFC > 0, ranked by |stat| descending, at most 250 genes;
    DOWN = genes with BH FDR < 0.05 and logFC < 0, same ranking and cap.
  Set names: MoTrPAC|species|tissue|layer|contrast|UP or DOWN (e.g. MoTrPAC|rat|SKM-GN|PROT|F_8w|UP).
  Phosphosite columns are collapsed to genes as in the store.

Loading
  GMT is the library format used by Enrichr and GSEA. Upload the .gmt file as a custom library in
  Enrichr-style tools, pass it to GSEA / fgsea / clusterProfiler::read.gmt / gseapy (gene_sets=path),
  or read it line by line (split on tabs; field 1 = name, field 2 = description, rest = genes).

Caveats
  Set membership depends on the FDR cutoff (0.05) and the cap (250); columns with few significant
  genes give small or empty sets, and a gene missing from a set may simply not be measured in that
  column. The full per-gene statistics are in the store (contrasts.parquet), not in the GMT.
  Built 2026-09-27T01:16:22 from store d615c2f0b9aec182.
