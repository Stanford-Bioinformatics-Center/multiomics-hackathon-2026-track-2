"""Metabolomics (METAB) as a third layer: store view, metabolite signatures and a 3-layer pathway panel.

Rows come from store/metab.parquet (store.build_metab; schema in store/SCHEMA.md "## METAB layer"). They have the
contrasts.parquet schema, but for layer == "METAB" `gene_symbol_human` holds the RefMet NAME of a metabolite.

Pathway panel: for each pathway in store/pathway_map.csv and each comparison column of one tissue, the cameraPR t
(legacy.camera_pr, inter-gene correlation 0.01) of the pathway's gene set (RNA, PROT columns) or metabolite set
(METAB columns) against every other feature measured in that column. t > 0 = the set moves up with exercise
relative to the rest of the column. The map is many-to-many: a pathway has a gene side and a metabolite side, each
tested within its own layer. Genes and metabolites are never matched one-to-one, and metabolite pools are not flux.
"""
import re
import sys

import numpy as np
import pandas as pd
from pandas.api.types import union_categoricals

from . import core, legacy, store
from .core import EXERCISE_KINDS

PATHWAYS = ["TCA cycle", "Glycolysis / gluconeogenesis", "Fatty-acid beta-oxidation", "BCAA catabolism",
            "Purine metabolism", "Pyrimidine metabolism"]
MIN_SET = 3                     # a layer x column cell needs >= 3 measured members, otherwise t = NaN
PANEL_LAYERS = ["RNA", "PROT", "METAB"]
METAB_ID_COLS = ("refmet_name", "hmdb", "kegg")
GENE_ID_COLS = ("gene_symbol", "uniprot", "ensembl", "rat_symbol")

PATHWAY_CAPTION = (
    "Cell = cameraPR t (limma cameraPR, inter-gene correlation 0.01) of the pathway's gene set (RNA, protein) or "
    "metabolite set (metabolites) against every other feature measured in that comparison column; t > 0 = the set "
    "moves up with exercise relative to the rest of the column, grey = fewer than 3 members measured. Row labels "
    "give the number of genes / metabolites measured (maximum over columns). Metabolite pools are not flux: a "
    "higher concentration can mean more production or less consumption. The pathway map is many-to-many "
    "(store/pathway_map.csv: gene side from msigdbr KEGG legacy sets, metabolite side from RefMet classes and "
    "curated names); genes and metabolites are never matched one-to-one."
)

# ------------------------------------------------------------------------------------------------ pathway map
# Gene side: one msigdbr gene set per pathway (store/pathway_genesets_msigdbr.csv, msigdbr 26.1.1, KEGG legacy).
GENE_SETS = {
    "TCA cycle": "KEGG_CITRATE_CYCLE_TCA_CYCLE",
    "Glycolysis / gluconeogenesis": "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
    "Fatty-acid beta-oxidation": "KEGG_FATTY_ACID_METABOLISM",
    "BCAA catabolism": "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
    "Purine metabolism": "KEGG_PURINE_METABOLISM",
    "Pyrimidine metabolism": "KEGG_PYRIMIDINE_METABOLISM",
}
# Metabolite side: RefMet sub_class rules and curated names (must exist in the MoTrPAC RefMet vocabulary).
# Every name is expanded to its notation variants (rat 'CAR(3:0)' == human 'CAR 3:0'; metab_features.name_key).
MET_CLASS = {
    "TCA cycle": ["TCA acids"],
    "Fatty-acid beta-oxidation": ["Acyl carnitines"],
}
MET_CLASS_NAME_RE = {"Fatty-acid beta-oxidation": r"^CAR[ (]"}   # acylcarnitines only (drops 3-Dehydroxycarnitine)
MET_CURATED = {
    "TCA cycle": ["Acetyl-CoA", "Succinyl-CoA"],
    "Glycolysis / gluconeogenesis": ["Glucose", "Hexose 6-phosphate", "Fructose 6-phosphate",
                                     "Fructose 1,6-bisphosphate", "Dihydroxyacetone phosphate",
                                     "Glyceraldehyde 3-phosphate", "Phosphoglyceric acid",
                                     "Phosphoenolpyruvic acid", "Pyruvic acid", "Lactic acid"],
    "Fatty-acid beta-oxidation": ["Carnitine"],
    "BCAA catabolism": ["Leucine", "Isoleucine", "Valine", "Leucine/Isoleucine", "Ketoleucine",
                        "3-Methyl-2-oxovaleric acid", "Ketoisovaleric acid", "Isovaleryl-CoA", "Propionyl-CoA",
                        "CAR(3:0)", "CAR(5:0)", "CAR(5:1)", "CAR(5:0(OH))"],
    "Purine metabolism": ["Adenine", "Guanine", "Hypoxanthine", "Xanthine", "Uric acid", "Allantoin",
                          "Adenosine", "Guanosine", "Inosine", "Xanthosine", "Deoxyadenosine", "Deoxyguanosine",
                          "AMP", "GMP", "IMP", "XMP", "Adenylsuccinic acid", "ADP", "GDP", "ATP", "GTP", "dATP",
                          "3',5' cyclic AMP", "Phosphoribosyl pyrophosphate"],
    "Pyrimidine metabolism": ["Uracil", "Cytosine", "Uridine", "Cytidine", "Thymidine", "Deoxycytidine",
                              "Deoxyuridine", "UMP", "CMP", "dUMP", "UDP", "CDP", "UTP", "CTP", "dCTP", "dTTP",
                              "Ureidopropionic acid", "beta-Alanine"],
}
CURATED_WHY = {
    "TCA cycle": "curated: TCA-cycle acyl-CoA (KEGG map00020 compound)",
    "Glycolysis / gluconeogenesis": "curated: glycolysis/gluconeogenesis intermediate (KEGG map00010 compound)",
    "Fatty-acid beta-oxidation": "curated: free carnitine (carnitine shuttle)",
    "BCAA catabolism": "curated: BCAA, branched-chain keto acid or BCAA-derived acyl-CoA/acylcarnitine "
                       "(KEGG map00280)",
    "Purine metabolism": "curated: purine base/nucleoside/nucleotide or degradation product (KEGG map00230)",
    "Pyrimidine metabolism": "curated: pyrimidine base/nucleoside/nucleotide or degradation product "
                             "(KEGG map00240)",
}


def _variants(feats):
    """RefMet name -> all RefMet names sharing its notation-insensitive key (incl. itself)."""
    by_key = feats.groupby("name_key").refmet_name.apply(list)
    return dict(zip(feats.refmet_name, feats.name_key.map(by_key)))


def build_pathway_map(write=True):
    """Write store/pathway_map.csv (pathway, side, member, rule, source, measured_in) and return it."""
    feats = store.load_metab_features()
    var = _variants(feats)
    sub = feats.set_index("refmet_name").sub_class
    kegg = feats.set_index("refmet_name").kegg
    msp = feats.set_index("refmet_name").species
    gs = pd.read_csv(store.STORE / "pathway_genesets_msigdbr.csv")
    rows = store.load(["layer", "gene_symbol_human"])
    glayers = rows.drop_duplicates().groupby("gene_symbol_human", observed=True).layer.apply(
        lambda s: "/".join(l for l in ["RNA", "PROT", "PHOSPHO", "PROT_OLINK"] if l in set(s)))
    out = []
    for pw in PATHWAYS:
        g = gs[gs.gs_name == GENE_SETS[pw]]
        assert len(g), GENE_SETS[pw]
        for r in g.drop_duplicates("gene_symbol").itertuples():
            out.append(dict(pathway=pw, side="gene", member=r.gene_symbol, rule=f"msigdbr gene set {r.gs_name}",
                            source=f"msigdbr {r.msigdbr_version} {r.collection} {r.gs_name} (KEGG {r.gs_exact_source})",
                            measured_in=glayers.get(r.gene_symbol, "")))
        mets = {}
        for sc in MET_CLASS.get(pw, []):
            names = feats.refmet_name[feats.sub_class == sc]
            if pw in MET_CLASS_NAME_RE:
                names = names[names.str.contains(MET_CLASS_NAME_RE[pw], regex=True)]
            rule = f"RefMet sub_class = {sc}" + (" (names CAR(...) / CAR ...)" if pw in MET_CLASS_NAME_RE else "")
            for n in names:
                mets.setdefault(n, rule)
        for n in MET_CURATED.get(pw, []):
            if n not in var:
                raise ValueError(f"curated {pw} metabolite {n!r} is not in the MoTrPAC RefMet vocabulary")
            for v in var[n]:
                mets.setdefault(v, CURATED_WHY[pw] + (f"; notation variant of {n}" if v != n else ""))
        for n, rule in mets.items():
            k = kegg.get(n)
            out.append(dict(pathway=pw, side="metabolite", member=n, rule=rule,
                            source=f"RefMet sub_class {sub.get(n)}" + (f"; KEGG {k}" if isinstance(k, str) and k
                                                                       else "") +
                                   " (store/raw_metab: rat_metab_feature_map_v1 / human_metab_feature_map_v1)",
                            measured_in=msp.get(n, "")))
    m = pd.DataFrame(out)
    m["_o"] = m.pathway.map({p: i for i, p in enumerate(PATHWAYS)})
    m = m.sort_values(["_o", "side", "member"]).drop(columns="_o").reset_index(drop=True)
    if write:
        m.to_csv(store.PATHWAY_MAP, index=False)
        store._CACHE.pop("pwmap", None)
    return m


def pathway_map():
    return store.load_pathway_map()


# ------------------------------------------------------------------------------------------------ store view
_CACHE = {}


def metab_store():
    """core.Store with the METAB rows appended and S.cols extended by the METAB column index (cached)."""
    if "S" in _CACHE:
        return _CACHE["S"]
    base = store.load()
    mrows, mcols = store.load_metab()
    merged = pd.concat([base, mrows], ignore_index=True)
    for c in base.columns:  # pd.concat drops to object when category sets differ; union them instead
        if isinstance(base[c].dtype, pd.CategoricalDtype):
            merged[c] = union_categoricals([base[c], mrows[c].astype("category")], ignore_order=True)
    # core.Store() reads its base rows from store.load(); hand it the merged frame through that cache for the
    # duration of the constructor so every derived array (gcode, stat, pos, ...) covers the METAB rows.
    key = None
    saved = store._CACHE.get(key)
    store._CACHE[key] = merged
    try:
        S = core.Store()
    finally:
        if saved is None:
            store._CACHE.pop(key, None)
        else:
            store._CACHE[key] = saved
    mc = mcols.set_index("column_id", drop=False)
    mc = mc[~mc.index.isin(S.cols.index)]
    S.cols = pd.concat([S.cols, mc[[c for c in S.cols.columns if c in mc.columns]]])
    _CACHE["S"] = S
    return S


def metab_columns(S, tissues=None):
    """METAB column_ids of exercise kinds (rat training vs sedentary, human exercise vs control), display order."""
    c = S.cols
    m = (c.layer == "METAB") & c.kind.isin(EXERCISE_KINDS)
    if tissues is not None:
        m &= c.tissue.isin([tissues] if isinstance(tissues, str) else tissues)
    return c[m].sort_values("col_order").column_id.tolist()


# ------------------------------------------------------------------------------------------------ signatures
HMDB_RE = re.compile(r"^HMDB(\d+)$", re.I)
KEGG_RE = re.compile(r"^(?:cpd:)?(C\d{5})$", re.I)


def _norm_hmdb(v):
    m = HMDB_RE.match(v.strip())
    return f"HMDB{int(m.group(1)):07d}" if m else None


def _norm_kegg(v):
    m = KEGG_RE.match(v.strip())
    return m.group(1).upper() if m else None


def is_metab_signature(path):
    """True if the CSV has a refmet_name / hmdb / kegg column and no gene id column."""
    try:
        cols = {c.strip().lower() for c in pd.read_csv(path, nrows=0).columns}
    except Exception:
        return False
    return bool(cols & set(METAB_ID_COLS)) and not cols & set(GENE_ID_COLS)


def map_metabolites(t):
    """Adds refmet (RefMet names joined by '; '), map_method, direction_parsed and status columns. A row maps to
    every notation variant of its RefMet name (rat 'CAR(2:0)' and human 'CAR 2:0'); HMDB / KEGG ids map to every
    RefMet name annotated with that id in metab_features."""
    from .signature import _dir
    feats = store.load_metab_features()
    var = _variants(feats)
    lower = {n.lower(): n for n in feats.refmet_name}
    key = feats.groupby("name_key").refmet_name.apply(list).to_dict()
    hm = feats.dropna(subset=["hmdb"]).groupby("hmdb").refmet_name.apply(list).to_dict()
    kg = feats.dropna(subset=["kegg"]).groupby("kegg").refmet_name.apply(list).to_dict()

    def expand(names):
        out = []
        for n in names:
            for v in var[n]:
                if v not in out:
                    out.append(v)
        return out

    def one(r):
        def val(c):
            v = r.get(c)
            return v.strip() if isinstance(v, str) and v.strip() else None
        n = val("refmet_name")
        if n:
            if n in var:
                return expand([n]), "RefMet name"
            if n.lower() in lower:
                return expand([lower[n.lower()]]), "RefMet name (case-normalised)"
            k = re.sub(r"[\s();]", "", n.lower())
            if k in key:
                return expand(key[k]), "RefMet name (notation-normalised)"
        h = val("hmdb")
        if h and _norm_hmdb(h) in hm:
            return expand(hm[_norm_hmdb(h)]), "HMDB id -> RefMet (metab_features)"
        kk = val("kegg")
        if kk and _norm_kegg(kk) in kg:
            return expand(kg[_norm_kegg(kk)]), "KEGG id -> RefMet (metab_features)"
        if n:
            return None, "unmapped: RefMet name not measured in MoTrPAC metabolomics"
        if h:
            return None, ("unmapped: HMDB id not in MoTrPAC metabolite annotation" if _norm_hmdb(h)
                          else "unmapped: not an HMDB id")
        if kk:
            return None, ("unmapped: KEGG id not in MoTrPAC metabolite annotation" if _norm_kegg(kk)
                          else "unmapped: not a KEGG compound id")
        return None, "unmapped: no refmet_name / hmdb / kegg value"

    got = [one(r) for r in t.to_dict("records")]
    t = t.copy()
    t["refmet"] = ["; ".join(g) if g else np.nan for g, _ in got]
    t["map_method"] = [m if g else np.nan for g, m in got]
    t["status"] = ["mapped" if g else m for g, m in got]
    t["direction_parsed"] = [_dir(v) for v in t.get("direction", pd.Series([None] * len(t)))]
    t.loc[(t.status == "mapped") & (t.direction_parsed == 0), "status"] = "dropped: direction not +1/-1"
    ok = t.status == "mapped"
    ex = t[ok].assign(name=t.refmet[ok].str.split("; ")).explode("name")
    dirs = ex.groupby("name").direction_parsed.nunique()
    conflict = set(dirs.index[dirs > 1])
    bad = ex[ex.name.isin(conflict)].index.unique()
    t.loc[bad, "status"] = "dropped: metabolite listed with conflicting directions"
    ok = t.status == "mapped"
    seen, dup = set(), []
    for i, s in t.refmet[ok].items():
        names = set(s.split("; "))
        if names <= seen:
            dup.append(i)
        seen |= names
    t.loc[dup, "status"] = "dropped: duplicate of an earlier row"
    return t


def load_metab_signature(path, name=None):
    """Read a metabolite signature CSV -> core.Sig whose `genes` are RefMet names (usable with core.score_columns
    on metab_store() and metab_columns()). Columns: refmet_name and/or hmdb / kegg, direction (+1/-1/up/down),
    optional group / weight / source."""
    from .signature import read_table
    t = read_table(path)
    if not set(METAB_ID_COLS) & set(t.columns):
        raise ValueError(f"{path}: needs a refmet_name, hmdb or kegg column")
    if "direction" not in t.columns:
        raise ValueError(f"{path}: needs a direction column (+1 up, -1 down)")
    if "group" not in t.columns:
        t["group"] = ""
    t["group"] = t.group.fillna("")
    t = map_metabolites(t)
    ok = t[t.status == "mapped"]
    genes, dirs, groups = [], [], {}
    for r in ok.itertuples():
        for n in r.refmet.split("; "):
            if n not in groups:
                genes.append(n)
                dirs.append(r.direction_parsed)
                groups[n] = r.group or "signature"
    dirs = np.asarray(dirs, int)
    return core.Sig(name=name or str(path).rsplit("/", 1)[-1].rsplit(".", 1)[0], table=t, genes=genes, dirs=dirs,
                    shown=list(genes), shown_dirs=dirs.copy(), groups=groups)


# ------------------------------------------------------------------------------------------------ pathway panel
def _panel_columns(S, tissue, species):
    ds = "human_acute" if species == "human" else "rat_train"
    kind = "exercise vs control" if species == "human" else "training vs sedentary"
    c = S.cols
    c = c[(c.dataset == ds) & (c.tissue == tissue) & (c.kind == kind) & c.layer.isin(PANEL_LAYERS)]
    per = c.groupby("contrast").layer.apply(set)
    keep = per.index[per.map(lambda s: set(PANEL_LAYERS) <= s)]
    return c[c.contrast.isin(keep)].sort_values("col_order")


def _time_label(c):
    return core.short_words(c)


def pathway_panel(S=None, tissue="SKM-GN", species="rat"):
    """Long table: pathway x (RNA|PROT|METAB) x column with cameraPR t / p of the pathway's gene set or metabolite
    set vs all other features measured in that column (t > 0 = set up with exercise). Only contrasts present in all
    three layers are used; t = NaN when fewer than MIN_SET members are measured."""
    S = metab_store() if S is None else S
    cols = _panel_columns(S, tissue, species)
    out_cols = ["pathway", "layer", "dataset", "tissue", "contrast", "time", "n_in_set", "n_universe", "t", "p",
                "column_id"]
    if not len(cols):
        return pd.DataFrame(columns=out_cols)
    pm = pathway_map()
    sets = {(pw, side): sorted(set(d.member)) for (pw, side), d in pm.groupby(["pathway", "side"])}
    rows = []
    for c in cols.itertuples():
        pos = S.col(c.column_id)
        stat = S.stat[pos]
        ok = np.isfinite(stat)
        g = S.gcode[pos][ok]
        stat = stat[ok]
        side = "metabolite" if c.layer == "METAB" else "gene"
        for pw in PATHWAYS:
            codes = S.codes(sets.get((pw, side), []))
            member = np.isin(g, codes[codes >= 0])
            m = int(member.sum())
            t = p = np.nan
            if m >= MIN_SET and m < len(stat) - 2:
                t, _, p = legacy.camera_pr(stat, member, core.INTER_GENE_COR)
            rows.append(dict(pathway=pw, layer=c.layer, dataset=c.dataset, tissue=c.tissue, contrast=c.contrast,
                             time=_time_label(c), n_in_set=m, n_universe=len(stat), t=float(t), p=float(p),
                             column_id=c.column_id))
    return pd.DataFrame(rows, columns=out_cols)


def pathway_heatmap(df, tissue):
    """Rows = pathway (measured genes / metabolites in the label), columns = layer x time; DIVERGING t, annotated."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    from .render import DIVERGING, INK, LAYER_COLORS, mpl_style
    layers_ = [l for l in PANEL_LAYERS if l in set(df.layer)]
    times = list(dict.fromkeys(df.time))  # df is in column display order
    pws = [p for p in PATHWAYS if p in set(df.pathway)]
    M = np.full((len(pws), len(layers_) * len(times)), np.nan)
    for r in df.itertuples():
        M[pws.index(r.pathway), layers_.index(r.layer) * len(times) + times.index(r.time)] = r.t
    nmax = df.groupby(["pathway", "layer"]).n_in_set.max()

    def lab(p):
        g = [f"{nmax.get((p, l), 0)} {core.LAYER_WORDS[l]}" for l in ["RNA", "PROT"] if l in layers_]
        mm = nmax.get((p, "METAB"), 0)
        return f"{p}\n{' / '.join(g)} genes · {mm} metabolites"
    vmax = max(3.0, float(np.nanmax(np.abs(M)))) if np.isfinite(M).any() else 3.0
    ncol = M.shape[1]
    with mpl_style():
        fig, ax = plt.subplots(figsize=(3.2 + 0.52 * ncol, 1.4 + 0.62 * len(pws)))
        im = ax.imshow(np.ma.masked_invalid(M), cmap=DIVERGING, norm=TwoSlopeNorm(0, -vmax, vmax), aspect="auto")
        for i in range(M.shape[0]):
            for j in range(ncol):
                v = M[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=10,
                            color="white" if abs(v) > 0.6 * vmax else INK)
        ax.set_yticks(range(len(pws)), [lab(p) for p in pws], fontsize=10)
        ax.set_xticks(range(ncol), times * len(layers_), rotation=90, fontsize=10)
        for k in range(1, len(layers_)):
            ax.axvline(k * len(times) - 0.5, color="white", lw=4)
        for k, l in enumerate(layers_):
            ax.text((k + 0.5) * len(times) - 0.5, -0.75, core.LAYER_WORDS[l], ha="center", va="bottom",
                    fontsize=11, color=LAYER_COLORS.get(l, INK), fontweight="bold")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
        cb.set_label("cameraPR t (+ = set up with exercise)", fontsize=10)
        cb.ax.tick_params(labelsize=10)
        sp = "Human" if (df.dataset == "human_acute").any() else "Rat"
        ax.set_title(f"{sp} {core.TISSUE_WORDS.get(tissue, tissue)}: pathway gene sets and metabolite sets",
                     fontsize=12, pad=26)
        fig.tight_layout()
    return fig


if __name__ == "__main__":  # python -m motrpac_probe.metab  -> rebuild store/pathway_map.csv
    m = build_pathway_map()
    print(m.groupby(["pathway", "side"]).size().to_string(), file=sys.stderr)
