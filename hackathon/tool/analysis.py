"""The core analysis in one readable file: does exercise push a disease signature back toward healthy, and which
omic layer (RNA or protein) shows it, and when?

    python analysis.py                                   # PAH skeletal-muscle signature (Malenfant 2015)
    python analysis.py examples/type2_diabetes_muscle_mootha2003.csv

Every step calls the motrpac_probe package; `mprobe run` runs the same steps and writes the full report.
"""
import sys

from motrpac_probe import core, layers, nulls, signature

SIG = sys.argv[1] if len(sys.argv) > 1 else "examples/pah_muscle_malenfant2015.csv"

# The comparisons that tell the story: one bout in humans (24 h after endurance or resistance exercise) versus
# weeks of training in rats (8 weeks), in RNA and in protein.
KEY = ["human_acute|VL|RNA|EE_vs_CON_24h", "human_acute|VL|PROT|EE_vs_CON_24h",
       "human_acute|VL|RNA|RE_vs_CON_24h", "human_acute|VL|PROT|RE_vs_CON_24h",
       "rat_train|SKM-GN|RNA|F_8w", "rat_train|SKM-GN|PROT|F_8w",
       "rat_train|SKM-GN|RNA|M_8w", "rat_train|SKM-GN|PROT|M_8w"]

# 1. The data. Every MoTrPAC comparison (human acute exercise, rat training; all tissues and layers) as one
#    gene-level table: logFC, test statistic and FDR per gene. Public summary statistics only.
S = core.get_store()

# 2. The question, written as a signature: genes up (+1) or down (-1) in the disease.
sig = signature.load(SIG)
print(f"Signature {sig.name}: {len(sig.genes)} genes counted "
      f"({(sig.table.status == 'mapped').sum()} of {len(sig.table)} rows mapped to MoTrPAC)\n")

# 3. Does exercise oppose it? Per comparison, count genes that exercise moves against their disease direction, and
#    test the whole set with cameraPR (limma's pre-ranked gene-set test). t > 0 means the set is opposed.
scores = core.score_columns(S, KEY, sig).set_index("column_id")

# 4. Is it specific? Repeat the set test on 1,000 random gene sets with the same pathway-class make-up
#    (mitochondrial / protein complex / secreted, and abundance). A percentile near 100 means the signature is
#    more opposed than its pathway class; a middling percentile means "any gene set like this would do".
null = nulls.run_nulls(S, KEY, sig, B=1000, kinds=("class",)).set_index("column_id")

print(f"{'comparison':<48}{'opposed':>9}{'cameraPR t':>12}{'class-null pct':>16}")
for cid in KEY:
    r, n = scores.loc[cid], null.loc[cid]
    print(f"{core.col_words(S.cols.loc[cid]):<48}{int(r.n_opposed):>4}/{int(r.n_measured):<4}{r.camera_t:>+12.2f}{n.pct_t:>16.0f}")

# 5. Which layer answered? For genes measured in both RNA and protein at the same time, compare the two layers and
#    the genome-wide RNA-protein agreement. After one bout RNA moves first; protein follows only with training.
pairs = layers.layer_pairs(S)
pairs = pairs[pairs.prot.isin(KEY) & ~pairs.cross]
summary, _ = layers.discordance(S, sig, pairs)
print()
for r in summary.itertuples():
    print(layers.sentence(r))

# 6. How to read this. 'Opposed' is a sign comparison between a disease-vs-healthy difference and an exercise
#    effect measured in other people or animals: a hypothesis, not evidence that exercise treats the disease.
print("\nGroup-level contrasts across genes, not per-person effects; see README 'Interpretation guardrails'.")
