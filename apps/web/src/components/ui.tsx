/* Small shared helpers for the explorer: number formatting, downloads, PNG export, rule switch. */

export type Significance = "raw" | "bh" | "bonferroni";
export const RULE_LABEL: Record<Significance, string> = { raw: "raw P < 0.05", bh: "BH q < 0.05", bonferroni: "Bonferroni < 0.05" };
const RULE_SHORT: Record<Significance, string> = { raw: "Raw P", bh: "BH q", bonferroni: "Bonferroni" };

export function fmt(value: number | null | undefined, digits = 3): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return "—";
  if (value === 0) return "0";
  if (Math.abs(value) < 0.001) return value.toExponential(1);
  return value.toFixed(digits).replace(/0+$/, "").replace(/\.$/, "");
}
export function signed(value: number | null | undefined, digits = 2): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return "—";
  return (value > 0 ? "+" : value < 0 ? "−" : "") + Math.abs(value).toFixed(digits);
}
export function download(name: string, contents: Blob): void {
  const link = document.createElement("a");
  const url = URL.createObjectURL(contents);
  link.href = url; link.download = name;
  document.body.append(link); link.click(); link.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function csvCell(value: unknown): string {
  const text = value === null || value === undefined ? "" : String(value);
  return '"' + text.replace(/"/g, '""') + '"';
}
export function csvDownload(name: string, rows: Array<Record<string, unknown>>): void {
  const columns = rows.length ? Object.keys(rows[0]) : ["empty"];
  const contents = [columns.map(csvCell).join(","), ...rows.map((row) => columns.map((c) => csvCell(row[c])).join(","))].join("\r\n");
  download(name, new Blob([contents], { type: "text/csv;charset=utf-8" }));
}
export async function svgPng(svg: SVGSVGElement | null, name: string): Promise<void> {
  if (!svg) return;
  const box = svg.viewBox.baseVal;
  const width = box && box.width ? box.width : svg.clientWidth || 900;
  const height = box && box.height ? box.height : svg.clientHeight || 400;
  const copy = svg.cloneNode(true) as SVGSVGElement;
  copy.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  copy.setAttribute("width", String(width));
  copy.setAttribute("height", String(height));
  copy.setAttribute("style", "font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif");
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)], { type: "image/svg+xml;charset=utf-8" }));
  const image = new Image();
  try {
    await new Promise<void>((resolve, reject) => { image.onload = () => resolve(); image.onerror = () => reject(new Error("render")); image.src = url; });
    const scale = 2, canvas = document.createElement("canvas");
    canvas.width = width * scale; canvas.height = height * scale;
    const context = canvas.getContext("2d");
    if (!context) return;
    context.scale(scale, scale);
    context.fillStyle = "#ffffff"; context.fillRect(0, 0, width, height);
    context.drawImage(image, 0, 0, width, height);
    const blob = await new Promise<Blob | null>((resolve) => canvas.toBlob(resolve, "image/png"));
    if (blob) download(name, blob);
  } finally { URL.revokeObjectURL(url); }
}
export function readFile(file: File | undefined, onText: (text: string) => void): void {
  if (file) file.text().then(onText).catch(() => undefined);
}

export function RuleSwitch({ rule, onChange }: { rule: Significance; onChange: (rule: Significance) => void }) {
  return <div className="lr-switch" role="group" aria-label="Significance rule">
    {(Object.keys(RULE_SHORT) as Significance[]).map((value) =>
      <button type="button" key={value} aria-pressed={rule === value} onClick={() => onChange(value)}
        title={value === "raw" ? "Unadjusted P < 0.05" : value === "bh" ? "Benjamini–Hochberg q < 0.05" : "Bonferroni-adjusted P < 0.05"}>
        {RULE_SHORT[value]}</button>)}
  </div>;
}

/** A pill that opens a small menu of options; unavailable options are shown but not selectable. */
export function Pick({ label, value, options, onChange }: {
  label: string; value: string; onChange: (value: string) => void;
  options: Array<{ value: string; label: string; disabled?: boolean; note?: string }>;
}) {
  const current = options.find((o) => o.value === value);
  return <details className={"lr-pill" + (value !== "all" ? " lr-pill--set" : "")}>
    <summary>{label}<b>{current?.label ?? value}</b><i aria-hidden="true">▾</i></summary>
    <div className="lr-menu">{options.map((o) =>
      <button type="button" key={o.value} aria-pressed={o.value === value} aria-disabled={o.disabled}
        className={"lr-menu-item" + (o.value === value ? " is-selected" : "") + (o.disabled ? " is-off" : "")}
        onClick={(event) => { if (o.disabled) return; onChange(o.value); event.currentTarget.closest("details")?.removeAttribute("open"); }}>
        {o.label}{o.note && <small>{o.note}</small>}
      </button>)}</div>
  </details>;
}
