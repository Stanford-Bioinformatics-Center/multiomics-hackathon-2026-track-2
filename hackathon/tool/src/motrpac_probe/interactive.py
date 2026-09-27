"""Vega-Lite spec for the interactive agreement grid (vendored Vega JS, inlined by render.py; no CDN)."""
import json

import numpy as np

from .core import LAYER_WORDS, TISSUE_WORDS, short_words


def grid_spec(S, L, cids, max_genes=80):
    c = S.cols.loc[list(cids)]
    info = {r.column_id: dict(block=f"{'Human' if r.dataset == 'human_acute' else 'Rat'} "
                                    f"{TISSUE_WORDS.get(r.tissue, r.tissue)}",
                              layer=LAYER_WORDS.get(r.layer, r.layer), when=short_words(r), order=int(i),
                              series=(r.group if r.dataset == "human_acute" else r.sex))
            for i, r in enumerate(c.itertuples())}
    genes = list(dict.fromkeys(L.gene))[:max_genes]
    vals = []
    for r in L[L.gene.isin(genes)].itertuples():
        d = info[r.column_id]
        f = lambda x, n=3: None if not np.isfinite(x) else round(float(x), n)
        vals.append(dict(gene=r.gene, direction="up in disease" if r.direction > 0 else "down in disease",
                         comparison=f"{d['block']} {d['layer']} {d['when']}", block=d["block"], layer=d["layer"],
                         when=d["when"], order=d["order"], series=d["series"], cell=f(r.cell, 2), stat=f(r.stat, 2),
                         logFC=f(r.logFC, 3), fdr_bh=f(r.fdr_bh, 4),
                         status="not measured" if not np.isfinite(r.stat) else
                         ("opposed" if r.agreement < 0 else "same direction")))
    blocks = sorted({v["block"] for v in vals})
    layers = sorted({v["layer"] for v in vals})
    spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {"values": vals},
        "params": [
            {"name": "tissue", "value": "all", "bind": {"input": "select", "options": ["all"] + blocks,
                                                        "name": "Tissue "}},
            {"name": "layerSel", "value": "all", "bind": {"input": "select", "options": ["all"] + layers,
                                                          "name": "Layer "}}],
        "transform": [{"filter": "(tissue == 'all' || datum.block == tissue) && "
                                 "(layerSel == 'all' || datum.layer == layerSel)"}],
        "vconcat": [
            {"width": {"step": 16}, "height": {"step": 16},
             "params": [{"name": "pick", "select": {"type": "point", "fields": ["gene"]},
                         "value": [{"gene": genes[0]}] if genes else None}],
             "mark": {"type": "rect", "stroke": "white", "strokeWidth": 1},
             "encoding": {
                 "x": {"field": "comparison", "type": "nominal", "sort": {"field": "order"}, "title": None,
                       "axis": {"labelAngle": -90, "labelFontSize": 11}},
                 "y": {"field": "gene", "type": "nominal", "sort": genes, "title": None},
                 "color": {"field": "cell", "type": "quantitative", "title": "agreement × |stat|",
                           "scale": {"domain": [-4, 4], "range": ["#1f4e79", "#f7f7f7", "#a51c30"], "clamp": True}},
                 "opacity": {"condition": {"param": "pick", "value": 1}, "value": 0.55},
                 "tooltip": [{"field": "gene"}, {"field": "direction"}, {"field": "comparison"},
                             {"field": "status"}, {"field": "logFC"}, {"field": "stat"}, {"field": "fdr_bh"}]}},
            {"width": 640, "height": 220,
             "transform": [{"filter": {"param": "pick"}}, {"filter": "isValid(datum.stat)"}],
             "mark": {"type": "line", "point": True},
             "encoding": {
                 "x": {"field": "when", "type": "ordinal", "sort": {"field": "order"}, "title": "comparison"},
                 "y": {"field": "stat", "type": "quantitative", "title": "gene stat (t or z)"},
                 "color": {"field": "layer", "type": "nominal", "title": "layer",
                           "scale": {"range": ["#2f5f8a", "#c26a2e", "#6b8e23", "#7a5195"]}},
                 "strokeDash": {"field": "series", "type": "nominal", "title": "group / sex"},
                 "detail": {"field": "block"},
                 "tooltip": [{"field": "gene"}, {"field": "comparison"}, {"field": "stat"}, {"field": "logFC"},
                             {"field": "fdr_bh"}]}}]}
    from .render import VEGA_CONFIG
    spec["config"] = VEGA_CONFIG
    return json.dumps(spec)
