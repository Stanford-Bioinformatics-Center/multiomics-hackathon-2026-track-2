"""Tests for motrpac_probe.render (HTML report rendering layer)."""
import json
import re
from html.parser import HTMLParser

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from motrpac_probe import render as R  # noqa: E402


# ---------------------------------------------------------------- formatting
@pytest.mark.parametrize(
    "p, expected",
    [
        (0.0001, "< 0.001"),
        (0.000999, "< 0.001"),
        (0.001, "0.001"),
        (0.003, "0.003"),
        (0.0034, "0.0034"),
        (0.034, "0.034"),
        (0.12, "0.12"),
        (0.5, "0.5"),
        (1.0, "1"),
        (np.float64(0.0456), "0.046"),
    ],
)
def test_fmt_p(p, expected):
    assert R.fmt_p(p) == expected


def test_fmt_p_na():
    assert R.fmt_p(float("nan")) == "NA"
    assert R.fmt_p(None) == "NA"


def test_fmt_num():
    assert R.fmt_num(1.234) == "1.23"
    assert R.fmt_num(1.234, digits=1) == "1.2"
    assert R.fmt_num(np.nan) == "not measured"
    assert R.fmt_num(None) == "not measured"
    assert R.fmt_num(pd.NA) == "not measured"
    assert R.fmt_num(12) == "12"
    assert R.fmt_num(0.00042) == "0.00042"


# ---------------------------------------------------------------- tables
def _df():
    return pd.DataFrame(
        {
            "Gene": ["A<b>", "B", "C"],
            "Effect": [1.5, np.nan, -0.25],
            "Adjusted p": [0.0001, 0.034, None],
        }
    )


def test_df_to_html_na_and_escape():
    html = R.df_to_html(_df(), table_id="t1", formats={"Adjusted p": R.fmt_p, "Effect": "{:.2f}"})
    assert html.count('<span class="na">not measured</span>') == 2
    assert "A&lt;b&gt;" in html and "A<b>" not in html
    assert "&lt; 0.001" in html and "0.034" in html and "1.50" in html
    assert 'class="data sortable"' in html and 'id="t1"' in html
    assert '<td class="num">' in html  # numeric columns right-aligned


def test_df_to_html_max_rows_and_row_class():
    df = pd.DataFrame({"x": range(10)})
    html = R.df_to_html(df, table_id="big", max_rows=3, row_class=lambda r: "flag" if r["x"] == 1 else None,
                        sortable=False)
    assert html.count("<tr") == 1 + 3
    assert "Showing 3 of 10 rows; the full table is in" in html
    assert "tables/big.csv" in html
    assert '<tr class="flag">' in html
    assert "sortable" not in html


# ---------------------------------------------------------------- figures
def _fig():
    with R.mpl_style():
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.imshow(np.arange(12).reshape(3, 4) - 6, cmap=R.DIVERGING)
        ax.plot([0, 1, 2], [0, 1, 0], color=R.LAYER_COLORS["RNA"])
        ax.set_title("Test figure")
    return fig


def test_fig_to_svg_viewbox_no_fixed_size(tmp_path):
    svg = R.save_fig(_fig(), tmp_path / "figures" / "f.png")
    assert (tmp_path / "figures" / "f.png").exists()
    assert svg.startswith("<svg")
    assert "<?xml" not in svg and "DOCTYPE" not in svg
    open_tag = re.match(r"<svg\b[^>]*>", svg).group(0)
    assert "viewBox=" in open_tag
    assert not re.search(r'\swidth="', open_tag)
    assert not re.search(r'\sheight="', open_tag)
    assert "<text" in svg  # text stays text
    assert "Test figure" in svg


def test_constants():
    assert set(R.LAYER_COLORS) == {"RNA", "PROT", "PHOSPHO", "METAB"}
    lo, hi = R.DIVERGING(0.0), R.DIVERGING(1.0)
    assert lo[2] > lo[0] and hi[0] > hi[2]  # blue negative, red positive


# ---------------------------------------------------------------- reports
def _context(tmp_path, *, interactive=False, fig_title="Opposition by tissue", how="Blue cells are opposed."):
    (tmp_path / "tables").mkdir(exist_ok=True)
    df = _df()
    df.to_csv(tmp_path / "tables" / "genes.csv", index=False)
    svg = R.save_fig(_fig(), tmp_path / "figures" / "heat.png")
    items = [
        {"kind": "figure", "title": fig_title, "how_to_read": how, "svg": svg,
         "png": "figures/heat.png", "id": "fig-heat"},
        {"kind": "table", "title": "Genes", "caption": "All genes.",
         "html": R.df_to_html(df, table_id="genes", formats={"Adjusted p": R.fmt_p}),
         "csv": "tables/genes.csv"},
        {"kind": "details", "summary": "Long table", "items": [
            {"kind": "table", "title": "Genes again", "caption": "",
             "html": R.df_to_html(df, table_id="genes2"), "csv": "tables/genes.csv"},
        ]},
    ]
    if interactive:
        spec = {"mark": "point", "data": {"values": [{"a": 1, "b": 2}]},
                "encoding": {"x": {"field": "a", "type": "quantitative"},
                             "y": {"field": "b", "type": "quantitative"},
                             "tooltip": [{"field": "a"}]}}
        items.append({"kind": "interactive", "title": "Scatter", "how_to_read": "Hover a point.",
                      "spec_json": json.dumps(spec), "id": "scatter", "fallback_svg": svg})
    return {
        "title": "Test report",
        "subtitle": "Minimal",
        "meta": [("Signature", "demo"), ("Genes", "3")],
        "headline": {"title": "Half the signature is opposed", "how_to_read": "Read it.",
                     "table_html": R.df_to_html(df.head(2), table_id="headline"),
                     "sentences": ["One sentence."]},
        "first_screen": ["One caveat."],
        "sections": [
            {"id": "main", "title": "Main results", "intro": ["Intro line."], "items": items},
            {"id": "hidden", "title": "Hidden section", "intro": [], "items": []},
        ],
        "toggles": {"hidden": False},
        "glossary": [("Opposed", "Opposite sign.")],
        "provenance_json": json.dumps({"a": 1}, indent=2),
    }


class _Refs(HTMLParser):
    def __init__(self):
        super().__init__()
        self.refs = []

    def handle_starttag(self, tag, attrs):
        for k, v in attrs:
            if k in ("src", "href"):
                self.refs.append((tag, v or ""))


def test_minimal_report_validates(tmp_path):
    out = R.render_report(_context(tmp_path), tmp_path / "report.html")
    assert R.validate_report(out) == []
    text = out.read_text(encoding="utf-8")
    # headline comes before the first section
    assert text.index("Half the signature is opposed") < text.index('id="main"')
    assert "Hidden section" not in text  # toggled off
    assert "<details>" in text and "How to read:" in text
    assert ":root" in text  # style.css inlined
    assert text.count("table.sortable") >= 1
    # self-contained: no http(s) script/link/img references
    p = _Refs()
    p.feed(text)
    for tag, url in p.refs:
        if tag in ("script", "link", "img"):
            assert not url.lower().startswith(("http", "//")), (tag, url)
    assert not re.search(r"<(script|link|img)[^>]+(src|href)=\"https?:", text)
    # no vega without interactive items
    assert "vegaEmbed" not in text and "vega-lite" not in text.lower()


def test_report_with_interactive_inlines_vega(tmp_path):
    out = R.render_report(_context(tmp_path, interactive=True), tmp_path / "report.html")
    text = out.read_text(encoding="utf-8")
    assert R.validate_report(out) == []
    if R._vendor_js() is not None:
        assert "vegaEmbed(" in text and "<noscript>" in text
        assert len(text) > 500_000
    assert not re.search(r"<script[^>]+src=", text)


def test_validate_catches_missing_caption(tmp_path):
    out = R.render_report(_context(tmp_path, how=""), tmp_path / "report.html")
    problems = R.validate_report(out)
    assert any("MISSING CAPTION" in p for p in problems)
    out = R.render_report(_context(tmp_path, fig_title=""), tmp_path / "report2.html")
    assert any("MISSING CAPTION" in p for p in R.validate_report(out))


def test_validate_catches_broken_image_external_script_empty_table(tmp_path):
    html = """<!DOCTYPE html><html><body>
    <img src="figures/nope.png" alt="x">
    <img src="figures/nope2.png">
    <script src="https://cdn.example.org/lib.js"></script>
    <link rel="stylesheet" href="http://example.org/s.css">
    <a href="https://doi.org/10.1000/xyz">citation</a>
    <a href="#top">top</a>
    <figure id="f1"><h3 class="fig-title">T</h3><svg viewBox="0 0 1 1"></svg></figure>
    <table><thead><tr><th>a</th></tr></thead><tbody></tbody></table>
    </body></html>"""
    path = tmp_path / "bad.html"
    path.write_text(html, encoding="utf-8")
    problems = R.validate_report(path)
    joined = "\n".join(problems)
    assert "figures/nope.png" in joined  # broken image ref
    assert "external <script>" in joined
    assert "external stylesheet" in joined
    assert "without alt" in joined
    assert "f1 has no how-to-read" in joined
    assert "empty table" in joined
    assert "doi.org" not in joined  # citation links allowed
    assert "#top" not in joined


def test_site_index(tmp_path):
    (tmp_path / "ex1").mkdir()
    (tmp_path / "ex1" / "report.html").write_text("<html></html>")
    rows = [
        {"name": "Example one", "href": "ex1/report.html", "description": "First", "n_genes": 12,
         "headline": {"RNA opposed": "61%", "p": "< 0.001"}},
        {"name": "Example two", "href": "ex1/report.html", "description": "Second", "n_genes": 5,
         "headline": {"RNA opposed": "40%"}},
    ]
    out = R.render_site_index(rows, tmp_path / "index.html", "Gallery", ["Intro."])
    text = out.read_text(encoding="utf-8")
    assert R.validate_report(out) == []
    assert "Example one" in text and "RNA opposed" in text and "&lt; 0.001" in text
    assert "not measured" in text  # missing headline value for row two
    rows[0]["href"] = "missing/report.html"
    out = R.render_site_index(rows, tmp_path / "index.html", "Gallery", [])
    assert any("missing/report.html" in p for p in R.validate_report(out))


def test_unknown_item_kind_raises(tmp_path):
    ctx = _context(tmp_path)
    ctx["sections"][0]["items"].append({"kind": "banner"})
    with pytest.raises(ValueError):
        R.render_report(ctx, tmp_path / "r.html")
