# Track 2: a reusable query engine and a testable discordance model

## Two separate claims

1. **Reusable query:** A processed disease-study signature can be mapped to compatible MoTrPAC feature, pathway, tissue, contrast, and time results. The output is a cross-study comparison with explicit coverage, mapping confidence, and incompatibility reasons. A shared sign does not mean the disease was treated by exercise.
2. **Track 2 model:** Within MoTrPAC human muscle, predict a held-out 24-hour protein effect from RNA and earlier phosphosite information, then ask whether those features improve out-of-fold prediction beyond simple baselines. Disease signatures are external examples and must not select the model or its thresholds.

This split matters because a large lookup table answers the catalog part of the brief but does not by itself answer its question, “Can this discordance be modeled or explained?” The [official human package overview](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/articles/package_overview.html) says it distributes “summary-level results only”; subject-level cross-omic covariance is therefore unavailable from that package.

## Input contract

The first app version accepts **processed feature-level results**. One row represents one feature and one named study contrast. Required fields are `study_id,species,tissue,layer,feature_id,id_namespace,log2_fc,statistic,p_value,q_value,contrast`; unavailable numeric values may be blank if the module does not require them. The meaning of a positive effect must be the named numerator minus denominator. Optional `biospecimen`, `gene_symbol`, `uniprot`, `refmet_id`, `refmet_name`, `source_feature_id`, and provenance/caveat fields stay attached to the result.

MoTrPAC reference results use the same fields plus `contrast_category` and `timepoint`. Keep a release version on every exported row and reject mixed releases in a run. The current local human reference is package 2.0.8, collection c2.0; the rat training reference uses its own source and intervention clock.

The Workbench importer is an **adapter with a required study-design manifest**. The [Workbench mwTab specification](https://www.metabolomicsworkbench.org/data/mwTab_specification.pdf) describes the subject identifier as optional and experimental factors as free name/value pairs. Thus neither pairing, comparator groups, nor confounder adjustment can be inferred safely from an arbitrary accession. [CFDE's NIH FAQ](https://commonfund.nih.gov/dataecosystem/faqs) says “The CFDE will not replicate or otherwise serve as a repository for these data sets”; use CFDE to discover studies and the Workbench download/API to obtain their data.

## Capability gates

| Input has | Safe analysis | Otherwise |
| --- | --- | --- |
| Stable gene, UniProt, or RefMet identifier | Exact or explicitly mapped feature context across measured MoTrPAC layers and times | Record unmatched or ambiguous identifier |
| Full measured feature universe and signed statistic | Rank association and a test using declared gene sets | Show queried hits only; do not infer enrichment |
| Named pathway with matching database/version | Pathway result comparison | Do not equate pathways by similar names |
| Explicit subject ID plus time labels | Paired change model, if design permits | Use an unpaired design only when justified; never invent pairs from sample numbering |
| Named metabolite identity | Direct metabolite lookup | Keep untargeted features as unmatched unless an audited mapping is supplied |
| Compatible tissue, species, and intervention contrast | Cross-study relation with caveat labels | Report the incompatibility; do not force a join |

The [Workbench REST specification](https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.2.pdf) offers study `factors`, `data`, and separate `untarg_data` endpoints. Retrieval can be generic; analysis eligibility is determined only after parsing and reviewing the study design.

## Model definition and validation

The primary outcome is the observed 24-hour MoTrPAC human-muscle total-protein log2 fold change for a uniquely mapped protein. The current catalog uses one preselected representative RNA and protein feature per gene or protein accession at the chosen target time; it retains RNA-only and protein-only entries and reports mapping exclusions. It distinguishes same-sign estimates, opposite-sign estimates, evidence-supported events, and indeterminate cases. A protein with q >= 0.05 is **not** called unchanged.

Compare three predictors on identical held-out proteins: zero change; RNA-only; RNA plus earlier RNA, earlier protein and parent-adjusted phosphosite summaries where available. The extra-feature model can only be compared on the same complete-case evaluation set; also report the coverage lost. Group all rows from one gene/protein into the same cross-validation fold. Report out-of-fold MAE and R²; if a binary discordance label can be defined with sufficient confidence, add calibration and prevalence-aware metrics. Repeat without each feature family to see whether timing or PTM information adds predictive value. Do not infer a causal mechanism from a feature's coefficient.

Two biological explanations are candidates: delayed protein response after RNA, and phosphorylation changes without a total-protein abundance shift. They are hypotheses. Human MoTrPAC muscle has three postexercise biopsy windows for these layers, which is too sparse for reliable gene-specific kinetic rates. The [MoTrPAC human muscle analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC13001352/) reports relatively few altered proteins at 24 hours and discusses temporally distinct RNA, protein, and phosphorylation responses. Rat training spans weeks and serves as a qualitative transfer check, not direct validation of an acute human predictor.

## What the current pilot actually found

The 24-hour EE-CON catalog has 15,411 entries, but 9,228 lack a paired total-protein measurement and 5,642 paired entries remain evidence-indeterminate under the declared thresholds. Only one pair is supported concordant, none supported opposite, and 285 show supported RNA response with a total-protein CI wholly within ±0.2 log2 FC. This large indeterminate class is a measurement/evidence result, not a biological statement that the layers agree.

For 5,586 unambiguous gene/protein pairs, gene-grouped out-of-fold R² for the later protein effect is 0.0101 using 24-hour RNA alone and 0.0116 after adding earlier RNA and site-minus-parent features. The richer model's MAE is essentially unchanged. EE→RE mode transfer has R² 0.0083. Thus the current code **does implement and evaluate a predictor**, but it does not yet explain much of the disagreement. Do not pitch the tiny R² increment as predictive success.

A stronger observed phenomenon is early phosphosite-versus-parent divergence: 735 phosphosite feature rows across 350 genes have an evidence-supported 15–45-minute site effect while their uniquely mapped parent-protein CI lies wholly within ±0.2. These are site rows, not 735 independent genes, and a site-minus-total-protein comparison is **not** measured phosphorylation occupancy. This supplies an explicit, measurable outcome for the next model.

For a stronger second model, restrict to uniquely mapped sites whose parent protein is within the equivalence band, then predict which sites show an evidence-supported early change using features available independently of that site's postexercise effect: baseline site/protein abundance, sequence motif, kinase or pathway annotation, and compartment. Split folds by **parent protein** so sites of one protein never appear in both training and test. Compare against prevalence and simple abundance baselines, and report precision–recall area, calibration, and coverage. Make any exercise-mode transfer an additional stress test; shared control subjects can limit independence. The 15–45-minute result has enough candidates to attempt this, whereas the later windows have only 13 and one such site rows under the same rules. This is a proposed next analysis, not a demonstrated prediction.

## What to show judges

1. One diagram of the study-design gate: upload → mapping/compatibility report → MoTrPAC context → model prediction.
2. One all-feature concordance/discordance catalog with the measured universe and ambiguous/unmatched counts.
3. One out-of-fold observed-versus-predicted protein plot and a baseline/ablation table. If the richer model does not improve prediction, report that result.
4. The early site-versus-parent plot and its mapped feature counts, with the occupancy caveat; it motivates the next predictive target.
5. One PAH example such as the nine resting-muscle proteins, labelled as external context. The earlier 24-hour RNA rise is exploratory; absent significant protein changes cannot establish that protein did not respond.

The central pitch should be narrow: **we built a reusable disease-to-MoTrPAC query and a measured discordance catalog, then asked whether timing and PTM features predict RNA–protein divergence. The first predictive test is weak, while early site-versus-parent divergence gives a stronger, testable next target.**
