"""Generalised MoTrPAC explorer: any molecule list in, every matching MoTrPAC result out.

Input lists carry no disease label. Each list is one kind:
    genes        gene symbols, UniProt, Ensembl or rat symbols -> RNA, protein, phospho, Olink layers
    metabolites  RefMet names (or HMDB / KEGG ids)              -> metabolomics layer
    pathways     MSigDB / MitoCarta / RefMet-class set names    -> MoTrPAC's precomputed CAMERA results
An optional `direction` (+1/-1, up/down) or numeric `score` column adds same/opposite-direction views and,
for scored gene lists, rank correlation.

Statistics are read from the motrpac_probe store (published MoTrPAC summary statistics). The only test run
here is the engine's cameraPR (core.score_column) for the whole input set in each comparison. One
multiple-testing family = every included MoTrPAC comparison of a layer (human and rat, all tissues, times,
contrasts); BH and Bonferroni are applied to it once, before any display filtering.
"""
from __future__ import annotations

import hashlib
import io
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from motrpac_probe import core, metab, signature
from motrpac_probe import store as _store_mod

API_ROOT = Path(__file__).resolve().parents[1]
REPO = API_ROOT.parents[1]
PATHWAY_TABLE = API_ROOT / "data" / "motrpac_camera_pathways.csv.gz"
MAX_ROWS_SHOWN = 40

# ---- examples: neutral names; the source is shown only as a citation ------------------------------
EXAMPLES = {
    "muscle_protein_9": dict(kind="genes", label="Muscle protein set · 9 proteins",
                             file=REPO / "MoTrPAC Hackathon/demo_inputs/3_Malenfant2015_Table2_muscle_down_proteins.xlsx",
                             source="Malenfant et al. 2015, J Mol Med, Table 2 (lower in the source group)"),
    "muscle_gene_87": dict(kind="genes", label="Muscle gene set · 87 genes",
                           file=REPO / "MoTrPAC Hackathon/demo_inputs/4_Mootha2003_T2D_muscle_OXPHOS_genes.xlsx",
                           source="Mootha et al. 2003, Nat Genet (MSigDB MOOTHA_VOXPHOS)"),
    "blood_gene_22": dict(kind="genes", label="Blood gene set · 22 genes",
                          file=REPO / "hackathon/tool/examples/pah_blood_cheadle2012_eds.csv",
                          source="Cheadle et al. 2012, PLoS ONE (erythroid subset of Table S2)"),
    "blood_ranked": dict(kind="genes", label="Blood ranked gene list · 19,593 genes",
                         file=REPO / "hackathon/tool/examples/pah_blood_gse33463_ranked.csv.gz",
                         source="GEO GSE33463, group comparison re-analysed with limma (moderated t)"),
    "blood_pathway_6": dict(kind="pathways", label="Blood pathway set · 8 pathways (6 BioCarta)",
                            file=REPO / "MoTrPAC Hackathon/demo_inputs/1_Cheadle2012_Table2_pathways.xlsx",
                            source="Cheadle et al. 2012, PLoS ONE, Table 2 selected pathways (lower in the source group)"),
    "plasma_metabolite_41": dict(kind="metabolites", label="Plasma metabolite set · 41 metabolites",
                                 file=REPO / "MoTrPAC Hackathon/demo_inputs/2_ST000763_plasma_metabolites.xlsx",
                                 source="Metabolomics Workbench ST000763, resting group difference (MoTrPAC Hackathon/Metabolomics)"),
    "tca_demo": dict(kind="metabolites", label="TCA intermediates · 10 metabolites (demo)",
                     file=REPO / "hackathon/tool/examples/tca_intermediates_demo.csv",
                     source="Synthetic demonstration list"),
}


def example_text(path) -> str:
    """An example file as CSV text; spreadsheets are read exactly as an upload would be."""
    path = Path(path)
    if path.suffix in (".xlsx", ".xls"):
        return read_upload(path.name, path.read_bytes())["text"]
    if path.suffix == ".gz":
        return pd.read_csv(path).to_csv(index=False)
    return path.read_text(encoding="utf-8")


def list_examples() -> list[dict]:
    return [{"id": k, "kind": v["kind"], "label": v["label"], "source": v["source"]}
            for k, v in EXAMPLES.items() if Path(v["file"]).is_file()]


# ---- MoTrPAC comparisons ---------------------------------------------------------------------------
LAYERS = {
    "transcriptomics": ("Transcriptomics", ["RNA"]),
    "proteomics": ("Proteomics", ["PROT", "PROT_OLINK"]),
    "phosphoproteomics": ("Phosphoproteomics", ["PHOSPHO"]),
    "metabolomics": ("Metabolomics", ["METAB"]),
}
TISSUE_LABEL = {
    "VL": "Skeletal muscle", "BLOOD": "Blood", "ADIPOSE": "Adipose", "SKM-GN": "Skeletal muscle (gastrocnemius)",
    "SKM-VL": "Skeletal muscle (vastus lateralis)", "HEART": "Heart", "LIVER": "Liver", "LUNG": "Lung",
    "KIDNEY": "Kidney", "WAT-SC": "Adipose (white, subcutaneous)", "BAT": "Adipose (brown)", "PLASMA": "Blood (plasma)",
    "ADRNL": "Adrenal", "COLON": "Colon", "CORTEX": "Brain cortex", "HIPPOC": "Hippocampus", "HYPOTH": "Hypothalamus",
    "OVARY": "Ovary", "TESTES": "Testes", "SMLINT": "Small intestine", "SPLEEN": "Spleen", "VENACV": "Vena cava",
}
HUMAN_TIME = {"20m-during": ("20 min (during)", 0), "40m-during": ("40 min (during)", 1), "10m": ("10 min after", 2),
              "15-45m": ("15–45 min after", 3), "3.5h": ("3.5–4 h after", 4), "24h": ("24 h after", 5)}
LAYER_LABEL = {"RNA": "RNA", "PROT": "Protein", "PROT_OLINK": "Plasma protein (Olink)", "PHOSPHO": "Phosphosite",
               "METAB": "Metabolite"}


def _category(c) -> str | None:
    if c.dataset == "rat_train":
        return "TRAIN-SED"
    k = str(c.contrast)
    if k.endswith("_pre"):
        return None
    for prefix, cat in (("EE_vs_CON_", "EE-CON"), ("RE_vs_CON_", "RE-CON"), ("EE_post_vs_pre_", "EE-EE"),
                        ("RE_post_vs_pre_", "RE-RE"), ("CON_post_vs_pre_", "CON-CON")):
        if k.startswith(prefix):
            return cat
    return None


EXERCISE = {"EE-CON": "endurance", "EE-EE": "endurance", "RE-CON": "resistance", "RE-RE": "resistance",
            "CON-CON": "control (no exercise)", "TRAIN-SED": "endurance training"}


@lru_cache(maxsize=1)
def _store():
    return metab.metab_store()


# ---- store-sourced study/dataset labels (no statistics computed here) -------------------------------
_DATASET_STUDY = {  # store dataset code -> (display phrase, store provenance version key)
    "rat_train": ("endurance training", "MotrpacRatTraining6moData"),
    "human_acute": ("acute exercise", "MotrpacHumanPreSuspensionAnalysis"),
}


def _dataset_for_species(species: str) -> str:
    """Store dataset code for a species, matching the store's own convention (human -> human_acute)."""
    return "human_acute" if species == "human" else "rat_train"


@lru_cache(maxsize=1)
def _provenance() -> dict:
    """Store provenance versions, e.g. {"MotrpacRatTraining6moData": "2.0.0", ...}; {} on any load failure."""
    try:
        return _store_mod.load_provenance()
    except Exception:
        return {}


def _study_label(dataset: str) -> str:
    """Human-readable study label sourced from store provenance; falls back to the dataset code (R4.7)."""
    entry = _DATASET_STUDY.get(dataset)
    if entry is None:
        return dataset  # unknown dataset -> the code itself
    phrase, version_key = entry
    version = _provenance().get(version_key)
    if not version:
        return dataset  # no store version -> the code itself, never fabricate
    return f"{phrase} ({version})"


def _columns(layer_codes: list[str]) -> pd.DataFrame:
    c = _store().cols.copy()
    c = c[c.layer.isin(layer_codes)]
    c["category"] = [_category(r) for r in c.itertuples()]
    c = c[c.category.notna()].copy()
    rat = c.dataset == "rat_train"
    c["time_key"] = np.where(rat, c.time.astype(str).str.replace("w", " wk"), c.time.astype(str))
    c["time_label"] = [t if r else HUMAN_TIME.get(t, (t, 9))[0] for t, r in zip(c.time_key, rat)]
    c["time_rank"] = [int(o) if r else HUMAN_TIME.get(str(t), ("", 9))[1] for t, o, r in zip(c.time, c.time_order, rat)]
    c["sex_label"] = np.where(rat, c.sex.map({"F": "female", "M": "male"}), "all")
    return c.sort_values(["species", "tissue", "layer", "category", "sex_label", "time_rank"])


# ---- input parsing ---------------------------------------------------------------------------------
ID_COLUMNS = {"genes": ["gene_symbol", "uniprot", "ensembl", "rat_symbol"],
              "metabolites": ["refmet_name", "hmdb", "kegg"],
              "pathways": ["pathway", "set", "set_name", "pathway_name", "paper_label", "name"]}
SOURCE_COLUMNS = ["source", "annotation", "database", "collection"]
SCORE_COLUMNS = ["score", "t", "stat", "logfc", "log2fc", "z"]


def detect_kind(columns: list[str]) -> str | None:
    """Guess the list kind from column names (used for uploaded spreadsheets)."""
    cols = " ".join(c.strip().lower() for c in columns)
    if any(k in cols for k in ("refmet", "metabolite", "hmdb", "kegg")):
        return "metabolites"
    if any(k in cols for k in ("pathway", "set_name", "biocarta", "gobp")) or re.search(r"set", cols):
        return "pathways"
    if any(k in cols for k in ("gene", "uniprot", "ensembl", "protein", "symbol")):
        return "genes"
    return None


def read_upload(filename: str, content: bytes) -> dict:
    """Uploaded CSV/TSV/TXT/XLSX -> {text (CSV), kind (guessed), n_rows}. The first sheet of a workbook is used."""
    name = filename.lower()
    if name.endswith((".xlsx", ".xlsm", ".xls")):
        t = pd.read_excel(io.BytesIO(content), sheet_name=0, dtype=str)
    else:
        t = pd.read_csv(io.BytesIO(content), sep=None, engine="python", dtype=str)
    t = t.dropna(how="all").dropna(axis=1, how="all")
    t.columns = [str(c).strip() for c in t.columns]
    if t.empty:
        raise ValueError("The file has no rows.")
    return {"text": t.to_csv(index=False), "kind": detect_kind(list(t.columns)), "n_rows": int(len(t)),
            "columns": list(t.columns)}


def _read(kind: str, text: str) -> pd.DataFrame:
    text = text.strip()
    if not text:
        raise ValueError("The list is empty.")
    first = text.splitlines()[0].lower()
    known = set(ID_COLUMNS[kind]) | {"direction", "score", "metabolite_name"} | set(SCORE_COLUMNS)
    header = any(h.strip() in known for h in re.split(r"[,\t]", first))
    if header:
        t = pd.read_csv(io.StringIO(text), sep=None, engine="python", dtype=str, keep_default_na=False, na_values=[""])
        t.columns = [c.strip().lower() for c in t.columns]
    else:  # bare list: one identifier per line, optional second field = direction or score
        rows = [re.split(r"[,\t;]", line.strip()) for line in text.splitlines() if line.strip()]
        t = pd.DataFrame({ID_COLUMNS[kind][0]: [r[0].strip() for r in rows]})
        if any(len(r) > 1 and r[1].strip() for r in rows):
            t["score"] = [r[1].strip() if len(r) > 1 else "" for r in rows]
    if kind == "metabolites" and "refmet_name" not in t.columns and "metabolite_name" in t.columns:
        t["refmet_name"] = t["metabolite_name"]
    if not set(ID_COLUMNS[kind]) & set(t.columns):
        raise ValueError(f"No identifier column found; use one of: {', '.join(ID_COLUMNS[kind])}.")
    score_col = next((c for c in SCORE_COLUMNS if c in t.columns), None)
    t["_score"] = pd.to_numeric(t[score_col], errors="coerce") if score_col else np.nan
    if "direction" in t.columns:
        t["_dir"] = [signature._dir(v) for v in t.direction]
    elif score_col:
        t["_dir"] = np.sign(t["_score"]).fillna(0).astype(int)
    else:
        t["_dir"] = 0
    return t


def _map(kind: str, t: pd.DataFrame):
    """(mapped table with name/dir/score, unmapped identifiers, directed, ranked, input rows mapped)."""
    directed = bool((t._dir != 0).any())
    ranked = bool(t._score.notna().sum() >= 50)
    work = t.copy()
    # The engine's mappers require a direction; undirected lists are mapped with a placeholder +1.
    work["direction"] = np.where(work._dir != 0, work._dir, 1) if directed else 1
    if kind == "genes":
        m = signature.map_genes(work, _store().genes)
        ok = m[m.status == "mapped"]
        rows_mapped = len(ok)
        out = pd.DataFrame({"name": ok.gene, "dir": ok._dir if directed else 0, "score": ok._score, "row": ok.index})
        unmapped = m.loc[m.status.str.startswith("unmapped"), [c for c in ID_COLUMNS["genes"] if c in m.columns]]
    else:
        m = metab.map_metabolites(work)
        rows_mapped = int((m.status == "mapped").sum())
        # One input row can match several MoTrPAC spellings of the same metabolite (human "CAR 3:0",
        # rat "CAR(3:0)"); `row` keeps them together as one molecule.
        ok = m[m.status == "mapped"].assign(name=lambda d: d.refmet.str.split("; "), row=lambda d: d.index).explode("name")
        out = pd.DataFrame({"name": ok.name, "dir": ok._dir if directed else 0, "score": ok._score, "row": ok.row})
        unmapped = m.loc[m.status.str.startswith("unmapped"), [c for c in ID_COLUMNS["metabolites"] if c in m.columns]]
    out = out.drop_duplicates("name").reset_index(drop=True)
    misses = [str(v) for v in unmapped.bfill(axis=1).iloc[:, 0].tolist()] if len(unmapped) else []
    return out, misses, directed, ranked, rows_mapped


# ---- per-layer analysis ----------------------------------------------------------------------------
def _family(p: np.ndarray):
    bh = np.full(len(p), np.nan)
    bonf = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        bh[ok] = false_discovery_control(p[ok], method="bh")
        bonf[ok] = np.minimum(p[ok] * ok.sum(), 1.0)
    return bh, bonf


def _layer(layer_id: str, mapped: pd.DataFrame, directed: bool, ranked: bool) -> dict:
    label, codes = LAYERS[layer_id]
    S = _store()
    cols = _columns(codes)
    names = set(S.genes)
    present = mapped[mapped.name.isin(names)]
    # Keep only molecules measured in at least one column of this layer; count where each spelling occurs.
    seen_all, seen_human = {}, {}
    for cid, species in zip(cols.column_id, cols.species):
        for name in S.genes[np.unique(S.gcode[S.col(cid)])]:
            seen_all[name] = seen_all.get(name, 0) + 1
            if species == "human":
                seen_human[name] = seen_human.get(name, 0) + 1
    present = present[present.name.isin(seen_all.keys())]
    if present.empty:
        return {"id": layer_id, "label": label, "available": False, "n_matched": 0,
                "reason": "None of the list's molecules is measured in this MoTrPAC layer."}
    # Undirected lists: every member gets d = -1, so the engine's signed cameraPR t reads "set shifts up".
    dirs = present.dir.to_numpy(int) if directed else -np.ones(len(present), int)
    sig = core.Sig(name="input", table=present, genes=present.name.tolist(), dirs=dirs)
    scores = pd.DataFrame([core.score_column(S, cid, sig) for cid in cols.column_id])
    if ranked:  # a whole ranked list is not a set; rank correlation below is its test
        for f in ("camera_t", "camera_p"):
            scores[f] = np.nan
    bh, bonf = _family(scores.camera_p.to_numpy(float))

    # One displayed molecule per input row; its label is the spelling used most in human data.
    molecules = []
    for row_id, group in present.groupby("row", sort=False):
        names = sorted(group.name, key=lambda n: (-seen_human.get(n, 0), -seen_all.get(n, 0), n))
        molecules.append({"name": names[0], "dir": int(group.dir.iloc[0]), "score": _f(group.score.iloc[0]),
                          "codes": [int(c) for c in S.codes(names) if c >= 0]})
    if ranked:
        molecules.sort(key=lambda m: -abs(m["score"] or 0))
    shown = molecules[:MAX_ROWS_SHOWN]
    score_map = dict(zip(present.name, present.score)) if ranked else {}
    columns, values = [], []
    for j, c in enumerate(cols.itertuples()):
        pos = S.col(c.column_id)
        g = S.gcode[pos]
        n_col = len(pos)
        s = scores.iloc[j]
        row = dict(id=c.column_id, species=c.species, tissue=c.tissue, tissue_label=TISSUE_LABEL.get(c.tissue, c.tissue),
                   layer=LAYER_LABEL.get(c.layer, c.layer), category=c.category, exercise=EXERCISE[c.category],
                   time=c.time_key, time_label=c.time_label, time_rank=int(c.time_rank), sex=c.sex_label,
                   dataset=c.dataset, study_label=_study_label(c.dataset),
                   n_tested=n_col, n_measured=int(s.n_measured),
                   set_t=_f(s.camera_t), set_p=_f(s.camera_p), set_bh=_f(bh[j]), set_bonferroni=_f(bonf[j]))
        if directed:
            row.update(n_opposite=int(s.n_opposed), n_same=int(s.n_same))
        else:
            row.update(n_up=int(s.n_opposed), n_down=int(s.n_same))
        if ranked:
            stat = pd.Series(S.stat[pos], index=S.genes[g])
            x = pd.Series(score_map).reindex(stat.index).to_numpy(float)
            rho, p, n = core.spearman(x, stat.to_numpy(float))
            row.update(rho=_f(rho), rho_p=_f(p), rho_n=n)
        columns.append(row)
        idx = {code: k for k, code in enumerate(g)}
        for i, mol in enumerate(shown):
            k = next((idx[code] for code in mol["codes"] if code in idx), None)
            if k is None:
                continue
            at = pos[k]
            p = float(S.rows.p.iat[at]) if "p" in S.rows else np.nan
            values.append([i, j, _f(S.lfc[at]), _f(p), _f(S.fdr[at]), _f(min(p * n_col, 1.0)) if np.isfinite(p) else None])
    if ranked:
        rho_p = np.array([c.get("rho_p") if c.get("rho_p") is not None else np.nan for c in columns], float)
        rbh, rbonf = _family(rho_p)
        for c, q, b in zip(columns, rbh, rbonf):
            c.update(rho_bh=_f(q), rho_bonferroni=_f(b))
    return {
        "id": layer_id, "label": label, "available": True, "n_matched": len(molecules),
        "molecules": [{k: v for k, v in m.items() if k != "codes"} for m in shown],
        "n_shown": int(len(shown)), "columns": columns, "values": values,
        "family": {"n_tests": int(np.isfinite(scores.camera_p.to_numpy(float)).sum()) if not ranked else len(columns),
                   "scope": f"every MoTrPAC {label.lower()} comparison (human and rat, all tissues, times and contrasts)",
                   "methods": ["BH", "Bonferroni"]},
        "molecule_rule": "Per-molecule BH = MoTrPAC BH within that comparison; Bonferroni = P × molecules tested in that comparison.",
    }


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


# ---- pathways ----------------------------------------------------------------------------------------
PATHWAY_TIME = {"during_20_min": ("20 min (during)", 0), "during_40_min": ("40 min (during)", 1),
                "post_10_min": ("10 min after", 2), "post_15_30_45_min": ("15–45 min after", 3),
                "post_3.5_4_hr": ("3.5–4 h after", 4), "post_24_hr": ("24 h after", 5)}
ASSAY_LABEL = {"transcript-rna-seq": "RNA", "prot-pr": "Protein", "prot-ph": "Phosphosite", "prot-ol": "Plasma protein (Olink)", "metab": "Metabolite"}


@lru_cache(maxsize=1)
def _pathways() -> pd.DataFrame:
    x = pd.read_csv(PATHWAY_TABLE).rename(columns={"z.std": "z"})
    fam = ["tissue", "assay", "contrast_short", "collection"]  # BH family verified to reproduce adj_p_value
    x["bonferroni"] = np.minimum(x.p_value * x.groupby(fam, observed=True).p_value.transform("size"), 1.0)
    parsed = x.contrast_short.str.extract(r"^(Endur|Resist|Control)\.(\S+) - (Endur|Resist|Control)\.(\S+)")
    cat = np.select([(parsed[0] == "Endur") & (parsed[2] == "Control"), (parsed[0] == "Resist") & (parsed[2] == "Control"),
                     (parsed[0] == "Endur") & (parsed[2] == "Endur"), (parsed[0] == "Resist") & (parsed[2] == "Resist"),
                     (parsed[0] == "Control") & (parsed[2] == "Control")],
                    ["EE-CON", "RE-CON", "EE-EE", "RE-RE", "CON-CON"], default="")
    x["category"], x["time"] = cat, parsed[1]
    x["key"] = x.set.str.upper()
    x["short_key"] = x.set_short.fillna("").str.upper().str.replace(r"[^A-Z0-9]", "", regex=True)
    return x[x.category != ""]


def _pathway_layer(t: pd.DataFrame, directed: bool) -> tuple[dict, list[str]]:
    if not PATHWAY_TABLE.is_file():
        return ({"id": "pathways", "label": "Pathways", "available": False, "n_matched": 0,
                 "reason": "MoTrPAC pathway table not built: run `Rscript apps/api/scripts_build_pathways.R` from the repo root."}, [])
    x = _pathways()
    id_col = next(c for c in ID_COLUMNS["pathways"] if c in t.columns)
    src_col = next((c for c in SOURCE_COLUMNS if c in t.columns), None)
    norm = lambda v: re.sub(r"[^A-Z0-9]", "", str(v).upper())
    sets = x.drop_duplicates("key")[["key", "database"]].copy()
    sets["full"] = sets.key.map(norm)
    sets["bare"] = [norm(k[len(d) + 1:]) if k.startswith(d + "_") else norm(k) for k, d in zip(sets.key, sets.database)]
    databases = set(sets.database)
    wanted, misses = [], []
    for i, (name, d) in enumerate(zip(t[id_col].astype(str), t._dir)):
        source = norm(t[src_col].iat[i]) if src_col and isinstance(t[src_col].iat[i], str) else None
        if source and source not in databases:
            misses.append(f"{name} ({t[src_col].iat[i]}: collection not in MoTrPAC)")
            continue
        pool = sets[sets.database == source] if source else sets
        key = norm(name)
        hit = pool[pool.full == key]
        if hit.empty:
            hit = pool[pool.bare == key]
        if hit.empty:
            misses.append(name)
            continue
        wanted.append((hit.key.iloc[0], int(d), name))
    if not wanted:
        return ({"id": "pathways", "label": "Pathways", "available": False, "n_matched": 0,
                 "reason": "No pathway name matched MoTrPAC's pathway results."}, misses)
    keys = [k for k, _, _ in wanted]
    sub = x[x.key.isin(keys)].copy()
    sub["time_label"] = sub.time.map(lambda v: PATHWAY_TIME.get(v, (v, 9))[0])
    sub["time_rank"] = sub.time.map(lambda v: PATHWAY_TIME.get(v, (v, 9))[1])
    cols = (sub[["tissue", "assay", "category", "time", "time_label", "time_rank"]].drop_duplicates()
            .sort_values(["tissue", "assay", "category", "time_rank"]).reset_index(drop=True))
    col_index = {tuple(r): j for j, r in enumerate(cols[["tissue", "assay", "category", "time"]].itertuples(index=False))}
    mol = {k: i for i, k in enumerate(keys)}
    values = [[mol[r.key], col_index[(r.tissue, r.assay, r.category, r.time)], _f(r.z), _f(r.p_value), _f(r.adj_p_value), _f(r.bonferroni)]
              for r in sub.itertuples()]
    set_short = sub.drop_duplicates("key").set_index("key").set_short.to_dict()
    columns = [dict(id=f"{c.tissue}|{c.assay}|{c.category}|{c.time}", species="human", tissue=c.tissue,
                    tissue_label={"muscle": "Skeletal muscle"}.get(c.tissue, c.tissue.capitalize()),
                    layer=ASSAY_LABEL.get(c.assay, c.assay), category=c.category, exercise=EXERCISE[c.category],
                    time=c.time, time_label=c.time_label, time_rank=int(c.time_rank), sex="all",
                    dataset=_dataset_for_species("human"), study_label=_study_label(_dataset_for_species("human"))) for c in cols.itertuples()]
    return ({"id": "pathways", "label": "Pathways", "available": True, "n_matched": len(keys),
             "molecules": [{"name": label if label.upper() != k else (set_short.get(k) or k), "dir": d if directed else 0, "score": None,
                            "set": k} for k, d, label in wanted],
             "n_shown": len(keys), "columns": columns, "values": values, "value_unit": "CAMERA z",
             "family": {"n_tests": None, "scope": "MoTrPAC's own pathway families (tissue × assay × contrast × collection)",
                        "methods": ["BH", "Bonferroni"]},
             "molecule_rule": "BH and Bonferroni within MoTrPAC's pathway family for that tissue, assay, contrast and collection."}, misses)


# ---- entry point -------------------------------------------------------------------------------------
_CACHE: dict[str, dict] = {}


def analyse(lists: list[dict]) -> dict:
    """lists: [{kind, name?, text? , example?}] -> inputs + one result per matched layer."""
    if not lists:
        raise ValueError("Add at least one list.")
    key = hashlib.sha256(repr([(l.get("kind"), l.get("example"), l.get("text")) for l in lists]).encode()).hexdigest()
    if key in _CACHE:
        return _CACHE[key]
    inputs, layers = [], []
    for i, item in enumerate(lists):
        kind = item.get("kind")
        if kind not in ID_COLUMNS:
            raise ValueError(f"Unknown list kind '{kind}'.")
        text, name, source = item.get("text"), item.get("name") or f"List {i + 1}", None
        if item.get("example"):
            ex = EXAMPLES.get(item["example"])
            if ex is None:
                raise ValueError(f"Unknown example '{item['example']}'.")
            text = example_text(ex["file"])
            name, source, kind = ex["label"], ex["source"], ex["kind"]
        t = _read(kind, text or "")
        if kind == "pathways":
            directed = bool((t._dir != 0).any())
            result, misses = _pathway_layer(t, directed)
            result["input"] = i
            layers.append(result)
            inputs.append(dict(id=i, name=name, kind=kind, source=source, n_rows=len(t), n_mapped=result["n_matched"],
                               unmapped=misses[:50], directed=directed, ranked=False))
            continue
        mapped, misses, directed, ranked, rows_mapped = _map(kind, t)
        inputs.append(dict(id=i, name=name, kind=kind, source=source, n_rows=len(t), n_mapped=rows_mapped, n_names=len(mapped),
                           unmapped=misses[:50], directed=directed, ranked=ranked))
        protein_list = kind == "genes" and any("uniprot" in c or "protein" in c for c in t.columns)
        order = (["metabolomics"] if kind == "metabolites" else
                 ["proteomics", "transcriptomics", "phosphoproteomics"] if protein_list else
                 ["transcriptomics", "proteomics", "phosphoproteomics"])
        for layer_id in order:
            result = _layer(layer_id, mapped, directed, ranked)
            result["input"] = i
            layers.append(result)
    layers.append({"id": "epigenomics", "label": "Epigenomics", "available": False, "n_matched": 0, "input": None,
                   "reason": "Methylation and chromatin results are not loaded into the explorer yet."})
    out = {"inputs": inputs, "layers": layers}
    _CACHE[key] = out
    return out
