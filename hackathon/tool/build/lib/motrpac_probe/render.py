"""HTML report rendering for motrpac_probe.

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
"""
from __future__ import annotations

import contextlib
import hashlib
import html as _html
import io
import json
import math
import re
from html.parser import HTMLParser
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping
from urllib.parse import unquote, urlsplit

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

__all__ = [
    "ACCENT", "INK", "MUTED", "DIVERGING", "LAYER_COLORS", "CYCLE", "MPL_RC", "VEGA_CONFIG",
    "fmt_p", "fmt_num", "df_to_html", "fig_to_svg", "save_fig", "mpl_style",
    "apply_mpl_style", "render_report", "render_site_index", "validate_report",
]

PKG_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = PKG_DIR / "templates"
VENDOR_DIR = PKG_DIR / "vendor"
VENDOR_FILES = ("vega.min.js", "vega-lite.min.js", "vega-embed.min.js")

# --------------------------------------------------------------------------
# Colours (figures). Keep in sync with --accent / --ink / --muted in style.css.
# --------------------------------------------------------------------------
ACCENT = "#2f5f8a"
INK = "#111111"
MUTED = "#6b6b6b"
GRID = "#d9d9d9"
#: Blue (negative) - near-white (0) - red (positive).
DIVERGING = LinearSegmentedColormap.from_list(
    "mprobe_diverging", ["#1f4e79", "#f7f7f7", "#a51c30"], N=256
)
DIVERGING.set_bad("#e6e6e6")
#: One restrained, colourblind-distinguishable colour per omic layer.
LAYER_COLORS = {
    "RNA": "#2f5f8a",
    "PROT": "#c26a2e",
    "PHOSPHO": "#6b8e23",
    "METAB": "#7a5195",
}
#: Default categorical cycle for figures (layer colours first, then greys).
CYCLE = [ACCENT, "#c26a2e", "#6b8e23", "#7a5195", "#8c8c8c", "#3b3b3b"]

#: Vega-Lite ``config`` block matching the report (pass as spec["config"]).
VEGA_CONFIG: dict[str, Any] = {
    "font": "Inter, Helvetica, Arial, sans-serif",
    "background": "white",
    "view": {"stroke": None},
    "axis": {"labelFontSize": 14, "titleFontSize": 15, "titleFontWeight": "normal",
             "labelColor": INK, "titleColor": INK, "domainColor": INK, "tickColor": INK,
             "grid": False},
    "legend": {"labelFontSize": 14, "titleFontSize": 14, "titleFontWeight": "normal"},
    "header": {"labelFontSize": 14, "titleFontSize": 15},
    "title": {"fontSize": 15, "fontWeight": "normal", "anchor": "start", "color": INK},
    "range": {"category": [ACCENT, "#c26a2e", "#6b8e23", "#7a5195", "#8c8c8c", "#3b3b3b"]},
    "mark": {"color": ACCENT},
}

MPL_RC: dict[str, Any] = {
    "figure.facecolor": "white",
    "figure.edgecolor": "white",
    "savefig.facecolor": "white",
    "savefig.edgecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": INK,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.titlesize": 12,
    "axes.titleweight": "normal",
    "axes.titlelocation": "left",
    "axes.labelsize": 11,
    "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
    "grid.color": GRID,
    "grid.linewidth": 0.5,
    "text.color": INK,
    "xtick.color": INK,
    "ytick.color": INK,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "legend.frameon": False,
    "font.size": 11,
    "font.family": "sans-serif",
    "font.sans-serif": ["Inter", "Helvetica", "Arial", "DejaVu Sans"],
    "svg.fonttype": "none",
    "svg.hashsalt": None,
    "image.cmap": "viridis",
    "lines.linewidth": 1.6,
    "patch.edgecolor": "white",
}


def mpl_style():
    """Return a context manager applying the report's matplotlib style.

    Usage::

        with mpl_style():
            fig, ax = plt.subplots(figsize=(8, 4))

    White background, black text, Inter/Helvetica/Arial with DejaVu Sans
    fallback, 11 pt base and 12 pt axes titles, no top/right spines, grid off,
    accent-first colour cycle, ``svg.fonttype='none'`` so SVG text stays text.
    See ``apply_mpl_style`` to set it globally instead.
    """
    return matplotlib.rc_context(MPL_RC)


def apply_mpl_style() -> None:
    """Apply the report's matplotlib style globally (``rcParams.update``)."""
    matplotlib.rcParams.update(MPL_RC)


# --------------------------------------------------------------------------
# Number formatting
# --------------------------------------------------------------------------
def _is_na(x: Any) -> bool:
    if x is None:
        return True
    try:
        res = pd.isna(x)
    except (TypeError, ValueError):
        return False
    return bool(res) if np.ndim(res) == 0 else False


def fmt_p(p: float) -> str:
    """Format a p value (or q value) for display.

    Returns ``"< 0.001"`` for p < 0.001, otherwise two significant digits
    (``0.0034`` never occurs; ``0.034``, ``0.12``, ``1``). NA input returns
    ``"NA"``, but callers should normally print their own NA text.
    """
    if _is_na(p):
        return "NA"
    p = float(p)
    if p < 0.001:
        return "< 0.001"
    s = f"{p:.2g}"
    if "e" in s:  # cannot happen for p >= 0.001, guard anyway
        s = f"{p:.4f}".rstrip("0").rstrip(".")
    return s


def fmt_num(x: Any, digits: int = 2) -> str:
    """Format a number with ``digits`` decimals; NA -> ``"not measured"``.

    Integers (and integer-valued numpy ints) print without decimals. Non-zero
    values that would round to zero print with ``digits`` significant digits
    instead (``0.00042`` rather than ``0.00``). Non-numeric input is returned
    as ``str(x)``.
    """
    if _is_na(x):
        return "not measured"
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if math.isinf(xf):
        return "inf" if xf > 0 else "-inf"
    if xf != 0 and abs(xf) < 0.5 * 10 ** (-digits):
        return f"{xf:.{max(digits, 1)}g}"
    return f"{xf:.{digits}f}"


def _default_cell(x: Any) -> str:
    if isinstance(x, (bool, np.bool_)):
        return "yes" if x else "no"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        xf = float(x)
        if xf.is_integer() and abs(xf) < 1e15:
            return str(int(xf))
        a = abs(xf)
        if a >= 100:
            return f"{xf:.0f}"
        if a >= 1:
            return f"{xf:.2f}"
        return f"{xf:.3g}"
    return str(x)


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------
def df_to_html(
    df: pd.DataFrame,
    *,
    table_id: str | None = None,
    caption: str | None = None,
    formats: Mapping[str, Callable[[Any], str] | str] | None = None,
    na: str = "not measured",
    sortable: bool = True,
    max_rows: int | None = None,
    number_cols_right: bool = True,
    row_class: Callable[[pd.Series], str | None] | None = None,
    html_cols: Iterable[str] = (),
) -> str:
    """Render a DataFrame as an HTML table string (wrapped in a scroll container).

    Parameters
    ----------
    df : DataFrame. Column names are used verbatim as headers (pass
        human-readable headers). The index is not rendered; reset it first
        if it carries information.
    table_id : ``id`` attribute; also names the full CSV in the truncation note.
    caption : optional ``<caption>`` text.
    formats : column -> callable (e.g. ``fmt_p``) or format string (``"{:.2f}"``).
        Formatters are not called on NA cells.
    na : text shown (in ``<span class="na">``) for NaN / None / pd.NA.
    sortable : add class ``sortable`` (the report's JS makes headers clickable).
    max_rows : show only the first ``max_rows`` rows plus a note linking to
        ``tables/<table_id>.csv`` (the caller must write that file).
    number_cols_right : right-align numeric columns.
    row_class : callable(row Series) -> CSS class or None (e.g. ``"flag"``).
    html_cols : columns whose values are trusted HTML (not escaped), e.g. links.

    All other cell values are HTML-escaped. An empty frame yields a table
    with a header and no body rows, which ``validate_report`` flags.
    """
    formats = dict(formats or {})
    html_cols = set(html_cols)
    n_total = len(df)
    shown = df if (max_rows is None or n_total <= max_rows) else df.iloc[:max_rows]

    numeric = {
        c: (pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c]))
        for c in df.columns
    }

    def fmt_cell(col: Any, val: Any) -> str:
        if _is_na(val):
            return f'<span class="na">{_html.escape(na)}</span>'
        f = formats.get(col)
        if f is None:
            s = _default_cell(val)
        elif callable(f):
            s = f(val)
        else:
            s = str(f).format(val)
        return s if col in html_cols else _html.escape(str(s))

    classes = ["data"] + (["sortable"] if sortable else [])
    id_attr = f' id="{_html.escape(str(table_id))}"' if table_id else ""
    out = [f'<div class="table-wrap"><table{id_attr} class="{" ".join(classes)}">']
    if caption:
        out.append(f"<caption>{_html.escape(str(caption))}</caption>")
    out.append("<thead><tr>")
    for c in df.columns:
        cls = ' class="num"' if (number_cols_right and numeric[c]) else ""
        out.append(f'<th scope="col"{cls}>{_html.escape(str(c))}</th>')
    out.append("</tr></thead><tbody>")
    cols = list(df.columns)
    for _, row in shown.iterrows():
        rc = row_class(row) if row_class else None
        out.append(f'<tr class="{_html.escape(rc)}">' if rc else "<tr>")
        for c in cols:
            cls = ' class="num"' if (number_cols_right and numeric[c]) else ""
            out.append(f"<td{cls}>{fmt_cell(c, row[c])}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    if len(shown) < n_total:
        if table_id:
            ref = f"tables/{table_id}.csv"
            where = f'<a href="{_html.escape(ref)}">{_html.escape(ref)}</a>'
        else:
            where = "the CSV file"
        out.append(
            f'<p class="note">Showing {len(shown)} of {n_total} rows; '
            f"the full table is in {where}</p>"
        )
    return "".join(out)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------
_RE_XML_HEAD = re.compile(r"<\?xml[^>]*\?>\s*", re.S)
_RE_DOCTYPE = re.compile(r"<!DOCTYPE[^>]*>\s*", re.S)
_RE_METADATA = re.compile(r"\s*<metadata>.*?</metadata>", re.S)
_RE_COMMENT = re.compile(r"<!--.*?-->", re.S)
_RE_SVG_OPEN = re.compile(r"<svg\b[^>]*>", re.S)


def fig_to_svg(fig) -> str:
    """Return a matplotlib figure as an inline-ready ``<svg>`` string.

    Text stays text (``svg.fonttype='none'``), the XML header, DOCTYPE,
    metadata and comments are stripped, fixed ``width``/``height`` attributes
    are removed (the ``viewBox`` keeps the native aspect ratio; CSS sets
    ``width: 100%``). A ``max-width`` equal to the native size (1 pt =
    4/3 px) is set inline so 11 pt text never renders larger than 11 pt.
    Element ids are prefixed with a content hash so several figures can be
    inlined in one page without id collisions. Does not close the figure.
    """
    buf = io.StringIO()
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig.savefig(buf, format="svg", bbox_inches="tight", metadata={"Date": None})
    svg = buf.getvalue()
    svg = _RE_XML_HEAD.sub("", svg)
    svg = _RE_DOCTYPE.sub("", svg)
    svg = _RE_METADATA.sub("", svg)
    svg = _RE_COMMENT.sub("", svg)

    m = _RE_SVG_OPEN.search(svg)
    if not m:
        raise ValueError("matplotlib did not produce an <svg> element")
    open_tag = m.group(0)
    vb = re.search(r'viewBox="([^"]+)"', open_tag)
    width_px = None
    if vb:
        parts = vb.group(1).split()
        if len(parts) == 4:
            width_px = float(parts[2]) * 4.0 / 3.0
    new_open = re.sub(r'\s(width|height)="[^"]*"', "", open_tag)
    style = f' style="max-width:{width_px:.0f}px"' if width_px else ""
    new_open = new_open.replace("<svg", f'<svg class="mpl-fig" role="img"{style}', 1)
    svg = svg[: m.start()] + new_open + svg[m.end():]

    # Scope matplotlib's global "*{...}" rule to this figure.
    svg = svg.replace("*{stroke-linejoin", ".mpl-fig *{stroke-linejoin")

    # Prefix ids to avoid collisions between figures in one document.
    pfx = "f" + hashlib.sha1(svg.encode()).hexdigest()[:8] + "-"
    svg = re.sub(r'\bid="([^"]+)"', lambda mm: f'id="{pfx}{mm.group(1)}"', svg)
    svg = re.sub(r"url\(#([^)]+)\)", lambda mm: f"url(#{pfx}{mm.group(1)})", svg)
    svg = re.sub(r'(xlink:href|href)="#([^"]+)"', lambda mm: f'{mm.group(1)}="#{pfx}{mm.group(2)}"', svg)
    return svg.strip()


def save_fig(fig, png_path: str | Path, dpi: int = 150) -> str:
    """Save ``fig`` as PNG (native aspect, ``bbox_inches='tight'``), return its inline SVG.

    Creates the parent directory if needed and closes the figure afterwards.
    """
    png_path = Path(png_path)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
        return fig_to_svg(fig)
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
_KINDS = {"figure", "table", "text", "interactive", "details"}


def _env():
    import jinja2

    return jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(TEMPLATE_DIR)),
        autoescape=jinja2.select_autoescape(["html", "j2", "html.j2"], default=True),
        trim_blocks=True,
        lstrip_blocks=True,
        undefined=jinja2.ChainableUndefined,
    )


def _walk_items(items: Iterable[Mapping]) -> Iterable[Mapping]:
    for it in items or []:
        yield it
        if it.get("kind") == "details":
            yield from _walk_items(it.get("items", []))


def _json_script_safe(s: str) -> str:
    """Make JSON text safe inside <script type="application/json"> ("<\\/" is a valid JSON escape)."""
    return s.replace("</", "<\\/").replace("<!--", "<\\u0021--")


def _js_script_safe(s: str) -> str:
    """Make JavaScript source safe inside a <script> element (escape only "</script")."""
    return re.sub(r"</(script)", r"<\\/\1", s, flags=re.I)


def _vendor_js() -> str | None:
    try:
        parts = [(VENDOR_DIR / f).read_text(encoding="utf-8") for f in VENDOR_FILES]
    except OSError:
        return None
    return "\n".join(_js_script_safe(p) for p in parts)


def _tool_version() -> str:
    try:
        from motrpac_probe import __version__

        return str(__version__)
    except Exception:  # pragma: no cover
        return "unknown"


def render_report(context: dict, out_html: Path) -> Path:
    """Render ``templates/report.html.j2`` with ``context`` and write ``out_html``.

    See the module docstring for the context contract. ``style.css`` and the
    table-sorting script are inlined; the vendored Vega, Vega-Lite and
    Vega-Embed builds are inlined only when an interactive item is present
    (if they are missing, interactive items show their fallback SVG).
    Raises ``ValueError`` for an unknown item kind. Returns the output path.
    """
    import markupsafe

    out_html = Path(out_html)
    ctx = dict(context)
    toggles = dict(ctx.get("toggles") or {})
    sections = [
        s for s in (ctx.get("sections") or []) if toggles.get(s.get("id"), True)
    ]
    counter = 0
    interactive_specs = []
    for s in sections:
        for it in _walk_items(s.get("items", [])):
            kind = it.get("kind")
            if kind not in _KINDS:
                raise ValueError(f"unknown report item kind: {kind!r} in section {s.get('id')!r}")
            if kind in ("figure", "interactive") and not it.get("id"):
                counter += 1
                it["id"] = f"{s.get('id', 'sec')}-fig{counter}"  # type: ignore[index]
            if kind == "interactive":
                interactive_specs.append(it)

    vendor = _vendor_js() if interactive_specs else None
    ctx["sections"] = sections
    ctx["show_glossary"] = bool(ctx.get("glossary")) and toggles.get("glossary", True)
    ctx["show_provenance"] = bool(ctx.get("provenance_json")) and toggles.get("provenance", True)
    ctx["toc"] = [(s["id"], s["title"]) for s in sections]
    if ctx["show_glossary"]:
        ctx["toc"].append(("glossary", "Glossary"))
    if ctx["show_provenance"]:
        ctx["toc"].append(("provenance", "Provenance"))
    ctx["css"] = markupsafe.Markup((TEMPLATE_DIR / "style.css").read_text(encoding="utf-8"))
    ctx["vega_js"] = markupsafe.Markup(vendor) if vendor else None
    ctx["has_interactive"] = bool(interactive_specs)
    ctx["script_safe"] = lambda s: markupsafe.Markup(_json_script_safe(str(s)))
    ctx.setdefault("tool_version", _tool_version())
    ctx.setdefault("headline", {})
    ctx.setdefault("meta", [])

    html_text = _env().get_template("report.html.j2").render(**ctx)
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html_text, encoding="utf-8")
    return out_html


def render_site_index(rows: list[dict], out_html: Path, title: str, intro: list[str]) -> Path:
    """Render the static gallery index (a plain table, one row per example).

    Each row: ``name`` (str), ``href`` (or ``report``; path to the example's report.html,
    relative to ``out_html``), ``description`` (str), ``n_genes`` (int or
    str), and ``headline``: an ordered mapping column -> value string. The
    headline columns are the union of all rows' keys, in first-seen order;
    missing values show "not measured". ``intro`` lines are plain text.
    """
    import markupsafe

    out_html = Path(out_html)
    cols: list[str] = []
    for r in rows:
        for k in (r.get("headline") or {}):
            if k not in cols:
                cols.append(k)
    html_text = _env().get_template("site_index.html.j2").render(
        title=title,
        intro=intro or [],
        rows=rows,
        headline_cols=cols,
        css=markupsafe.Markup((TEMPLATE_DIR / "style.css").read_text(encoding="utf-8")),
        tool_version=_tool_version(),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html_text, encoding="utf-8")
    return out_html


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------
_EXTERNAL_SRC_TAGS = {"script", "img", "iframe", "embed", "source", "video", "audio", "object", "input"}


def _is_external(url: str) -> bool:
    u = url.strip().lower()
    return u.startswith(("http://", "https://", "//"))


class _ReportParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.refs: list[tuple[str, str, str]] = []  # (tag, attr, value)
        self.problems: list[str] = []
        self.figures: list[dict] = []  # closed figures
        self._fig_stack: list[dict] = []
        self._how_depth: list[tuple[str, dict]] = []  # (tag, figure record)
        self._tables: list[dict] = []
        self.n_tables = 0
        self.empty_tables: list[str] = []
        self.n_missing = 0
        self.ids: dict[str, int] = {}
        self._svg_depth = 0  # inside <svg> or <noscript>: ids there are not checked
        self._raw_depth = 0  # inside <style>/<script>

    def _start(self, tag: str, attrs: list[tuple[str, str | None]], selfclosing: bool) -> None:
        a = {k.lower(): (v or "") for k, v in attrs}
        if tag in ("style", "script") and not selfclosing:
            self._raw_depth += 1
        if a.get("id") and not self._svg_depth:
            self.ids[a["id"]] = self.ids.get(a["id"], 0) + 1
        if tag in ("svg", "noscript") and not selfclosing:
            self._svg_depth += 1
        for attr in ("src", "href", "xlink:href", "data"):
            if attr in a and a[attr]:
                self.refs.append((tag, attr, a[attr]))
        if tag == "link" and a.get("href") and _is_external(a["href"]):
            self.problems.append(f"external stylesheet/link reference: {a['href']}")
        if tag in _EXTERNAL_SRC_TAGS and a.get("src") and _is_external(a["src"]):
            self.problems.append(f"external <{tag}> reference: {a['src']}")
        if tag == "img" and not a.get("alt", "").strip():
            self.problems.append(f"<img> without alt text: {a.get('src', '')}")
        classes = a.get("class", "").split()
        if tag == "figure" and not selfclosing:
            rec = {"id": a.get("id", f"(line {self.getpos()[0]})"), "how": False}
            self._fig_stack.append(rec)
        if "how" in classes and self._fig_stack and not selfclosing:
            self._how_depth.append((tag, self._fig_stack[-1]))
        if tag == "table" and not selfclosing:
            self._tables.append({"id": a.get("id") or f"(line {self.getpos()[0]})", "cells": 0})
        if tag == "td" and self._tables:
            self._tables[-1]["cells"] += 1

    def handle_starttag(self, tag, attrs):
        self._start(tag, attrs, False)

    def handle_startendtag(self, tag, attrs):
        self._start(tag, attrs, True)

    def handle_endtag(self, tag):
        if tag in ("style", "script") and self._raw_depth:
            self._raw_depth -= 1
        if tag in ("svg", "noscript") and self._svg_depth:
            self._svg_depth -= 1
        if self._how_depth and self._how_depth[-1][0] == tag:
            self._how_depth.pop()
        if tag == "figure" and self._fig_stack:
            self.figures.append(self._fig_stack.pop())
        if tag == "table" and self._tables:
            t = self._tables.pop()
            self.n_tables += 1
            if t["cells"] == 0:
                self.empty_tables.append(t["id"])

    def handle_data(self, data):
        if self._raw_depth:
            return
        self.n_missing += data.count("MISSING CAPTION")
        if self._how_depth and data.strip():
            _, rec = self._how_depth[-1]
            text = data.strip()
            if text.lower().rstrip(":") not in ("how to read",):
                rec["how"] = True


def validate_report(html_path) -> list[str]:
    """Check a rendered report (or site index) and return a list of problems.

    Checks: relative ``src``/``href`` targets exist next to the HTML file
    (``#anchors``, ``data:``, ``mailto:`` and http(s) links are skipped);
    no "MISSING CAPTION" marker; every ``<figure>`` has a non-empty
    how-to-read line (``class="how"``); every ``<img>`` has ``alt`` text; no
    external http(s) ``<script>``/``<link>``/``<img>`` (or other embedded
    resource) references (plain ``<a href="https://...">`` is allowed); no
    tables without body cells; no duplicate HTML ``id`` values (ids inside
    inline SVG are ignored). An empty list means the report passed.
    """
    html_path = Path(html_path)
    text = html_path.read_text(encoding="utf-8")
    p = _ReportParser()
    p.feed(text)
    p.close()
    problems = list(p.problems)

    n_missing = p.n_missing
    if n_missing:
        problems.append(f"MISSING CAPTION appears {n_missing} time(s)")

    base = html_path.parent
    seen = set()
    for tag, attr, val in p.refs:
        v = val.strip()
        if not v or v.startswith("#"):
            continue
        low = v.lower()
        if low.startswith(("http://", "https://", "//", "data:", "mailto:", "javascript:", "tel:")):
            continue
        if re.match(r"^[a-z][a-z0-9+.-]*:", low):  # other schemes
            continue
        path = unquote(urlsplit(v).path)
        if not path or (tag, path) in seen:
            continue
        seen.add((tag, path))
        if not (base / path).exists():
            problems.append(f"broken local reference <{tag} {attr}=\"{v}\">")

    for fig in p.figures:
        if not fig["how"]:
            problems.append(f"figure {fig['id']} has no how-to-read line")
    for t in p.empty_tables:
        problems.append(f"empty table {t}")
    for i, n in p.ids.items():
        if n > 1:
            problems.append(f"duplicate id \"{i}\" ({n} elements)")
    return problems
