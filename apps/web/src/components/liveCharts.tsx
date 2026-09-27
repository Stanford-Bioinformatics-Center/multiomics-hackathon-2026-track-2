import type { RefObject } from "react";

/* Presentation charts for section 05. Every chart is one SVG with a viewBox, so it scales to its
   card and exports to PNG as drawn. Diverging colour: burgundy = higher, slate blue = lower,
   warm grey = near zero. Significance is carried by a white dot, never by colour alone. */

const WARM = [140, 39, 64];   // #8c2740
const COOL = [58, 98, 140];   // #3a628c
const MID = [240, 237, 232];  // #f0ede8
const INK = "#33312e", SOFT = "#6d6862", FAINT = "#9a938b", RULE = "#e2ddd6";
export const COLOURS = { warm: "#8c2740", cool: "#3a628c", mid: "#f0ede8", rule: RULE, ink: INK, soft: SOFT };

function mix(a: number[], b: number[], t: number): string {
  const c = a.map((v, i) => Math.round(v + (b[i] - v) * t));
  return `rgb(${c[0]}, ${c[1]}, ${c[2]})`;
}
export function diverging(value: number | null | undefined, scale: number): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return "#faf9f7";
  const t = Math.max(-1, Math.min(1, value / (scale || 1)));
  return t >= 0 ? mix(MID, WARM, Math.sqrt(t)) : mix(MID, COOL, Math.sqrt(-t));
}
function short(value: number): string {
  if (!Number.isFinite(value)) return "—";
  if (Math.abs(value) >= 10) return value.toFixed(0);
  if (Math.abs(value) >= 1) return value.toFixed(1);
  return value.toFixed(2).replace(/^(-?)0\./, "$1.");
}

export type Cell = { value: number | null; pass: boolean | null; title: string } | null;
export type Col = { key: string; label: string };

function Legend({ x, y, scale, labels, dotKey = true }: { x: number; y: number; scale: number; labels: [string, string]; dotKey?: boolean }) {
  const steps = [-1, -0.5, 0, 0.5, 1];
  return <g transform={`translate(${x}, ${y})`}>
    <text x="0" y="12" fontSize="13" fill={SOFT} textAnchor="end">{labels[0]}</text>
    {steps.map((s, i) => <rect key={s} x={10 + i * 30} y="0" width="28" height="16" rx="4" fill={diverging(s * scale, scale)} />)}
    <text x={10 + steps.length * 30 + 6} y="12" fontSize="13" fill={SOFT}>{labels[1]}</text>
    {dotKey && <circle cx={10 + steps.length * 30 + 6 + labels[1].length * 7.2 + 26} cy="8" r="5" fill="#fff" stroke={FAINT} />}
    {dotKey && <text x={10 + steps.length * 30 + 6 + labels[1].length * 7.2 + 36} y="12" fontSize="13" fill={SOFT}>significant · faded: not significant</text>}
  </g>;
}

export function HeatGrid({ rows, groups, cell, scale, minScale = 0, labels, svgRef, rowLabelWidth = 150, onCell }: {
  rows: string[]; groups: Array<{ title?: string; cols: Col[] }>; cell: (row: string, col: string) => Cell;
  scale?: number; minScale?: number; labels: [string, string]; svgRef?: RefObject<SVGSVGElement | null>; rowLabelWidth?: number;
  onCell?: (row: string, col: string) => void;
}) {
  const nCols = groups.reduce((n, g) => n + g.cols.length, 0);
  const gap = 4, groupGap = 22, cellH = rows.length > 6 ? 34 : 44;
  const cellW = Math.min(118, (960 - rowLabelWidth - 24 - groupGap * (groups.length - 1)) / nCols - gap);
  const hasTitles = groups.some((g) => g.title);
  const top = hasTitles ? 58 : 34;
  const width = rowLabelWidth + nCols * (cellW + gap) + groupGap * (groups.length - 1) + 8;
  const height = top + rows.length * (cellH + gap) + 58;
  const values = rows.flatMap((r) => groups.flatMap((g) => g.cols.map((c) => cell(r, c.key)?.value ?? 0)));
  const s = scale ?? Math.max(minScale, 1e-9, ...values.map((v) => Math.abs(v ?? 0)));
  let x0 = rowLabelWidth;
  const placed = groups.map((g) => { const start = x0; x0 += g.cols.length * (cellW + gap) + groupGap; return { ...g, start }; });
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg"
    aria-label={`Heat map, ${rows.length} rows by ${nCols} columns`} style={{ maxWidth: width }}>
    {placed.map((g) => <g key={g.start}>
      {g.title && <text x={g.start} y="18" fontSize="14" fontWeight="600" fill={INK}>{g.title}</text>}
      {g.cols.map((c, i) => <text key={c.key} x={g.start + i * (cellW + gap) + cellW / 2} y={top - 12}
        fontSize="12.5" fill={SOFT} textAnchor="middle">{c.label}</text>)}
      {rows.map((r, j) => g.cols.map((c, i) => {
        const v = cell(r, c.key);
        const x = g.start + i * (cellW + gap), y = top + j * (cellH + gap);
        return <g key={r + c.key} onClick={onCell && v ? () => onCell(r, c.key) : undefined} style={onCell && v ? { cursor: "pointer" } : undefined}>
          <title>{v?.title ?? `${r} · ${c.label}: not measured`}</title>
          <rect x={x} y={y} width={cellW} height={cellH} rx="7" fill={v ? diverging(v.value, s) : "#faf9f7"}
            opacity={v?.pass === false ? 0.38 : 1} stroke={v ? "none" : RULE} strokeDasharray={v ? undefined : "3 3"} />
          {v?.pass && <circle cx={x + cellW / 2} cy={y + cellH / 2} r="5" fill="#fff" />}
        </g>;
      }))}
    </g>)}
    {rows.map((r, j) => <text key={r} x={rowLabelWidth - 14} y={top + j * (cellH + gap) + cellH / 2 + 5}
      fontSize="14" fill={INK} textAnchor="end">{r}</text>)}
    <Legend x={rowLabelWidth + 120} y={height - 28} scale={s} labels={labels} />
  </svg>;
}

export function Bars({ items, unit, svgRef, labels }: {
  items: Array<{ label: string; value: number | null; pass: boolean | null; title: string }>;
  unit: string; svgRef?: RefObject<SVGSVGElement | null>; labels?: [string, string];
}) {
  const width = 900, height = 340, top = 44, bottom = 262, left = 70;
  const max = Math.max(1e-9, ...items.map((i) => Math.abs(i.value ?? 0)));
  const hasNeg = items.some((i) => (i.value ?? 0) < 0), hasPos = items.some((i) => (i.value ?? 0) > 0);
  const zero = hasNeg && hasPos ? (top + bottom) / 2 : hasNeg ? top : bottom;
  const span = hasNeg && hasPos ? (bottom - top) / 2 : bottom - top;
  const slot = (width - left - 20) / items.length, barW = Math.min(64, slot * 0.46);
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg" aria-label={unit} style={{ maxWidth: 820 }}>
    <text x={left} y="18" fontSize="13" fill={SOFT}>{unit}</text>
    <line x1={left} x2={width - 20} y1={zero} y2={zero} stroke="#bdb6ae" />
    {items.map((item, i) => {
      const v = item.value ?? 0, h = Math.abs(v) / max * span, cx = left + slot * i + slot / 2;
      const y = v >= 0 ? zero - h : zero;
      const colour = v >= 0 ? COLOURS.warm : COLOURS.cool;
      return <g key={item.label}>
        <title>{item.title}</title>
        <rect x={cx - barW / 2} y={y} width={barW} height={Math.max(h, 1.5)} rx="6" fill={colour} opacity={item.pass === false ? 0.3 : 1} />
        <text x={cx} y={v >= 0 ? y - 10 : y + h + 22} fontSize="15" fontWeight="600" fill={INK} textAnchor="middle">{short(v)}</text>
        <text x={cx} y={height - 36} fontSize="13" fill={SOFT} textAnchor="middle">{item.label}</text>
      </g>;
    })}
    <text x={left} y={height - 8} fontSize="12.5" fill={FAINT}>
      {labels ? `${labels[1]} above the line · ${labels[0]} below · ` : ""}faded: not significant under the selected rule
    </text>
  </svg>;
}

export function PairBars({ items, names, svgRef }: {
  items: Array<{ label: string; a: number; b: number }>; names: [string, string]; svgRef?: RefObject<SVGSVGElement | null>;
}) {
  const width = 900, height = 360, top = 50, bottom = 290, left = 40;
  const max = Math.max(1, ...items.flatMap((i) => [i.a, i.b]));
  const slot = (width - left - 20) / items.length, barW = Math.min(34, slot * 0.28);
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg" aria-label="Paired counts by time" style={{ maxWidth: 820 }}>
    <circle cx={left + 6} cy="16" r="6" fill={COLOURS.warm} /><text x={left + 18} y="21" fontSize="13.5" fill={INK}>{names[0]}</text>
    <circle cx={left + 206} cy="16" r="6" fill={COLOURS.cool} /><text x={left + 218} y="21" fontSize="13.5" fill={INK}>{names[1]}</text>
    <line x1={left} x2={width - 20} y1={bottom} y2={bottom} stroke="#bdb6ae" />
    {items.map((item, i) => {
      const cx = left + slot * i + slot / 2;
      return <g key={item.label}>
        {[item.a, item.b].map((v, k) => {
          const h = v / max * (bottom - top), x = cx + (k ? 3 : -barW - 3);
          return <g key={k}>
            <rect x={x} y={bottom - h} width={barW} height={Math.max(h, 1.5)} rx="5" fill={k ? COLOURS.cool : COLOURS.warm} />
            <text x={x + barW / 2} y={bottom - h - 8} fontSize="14" fontWeight="600" fill={INK} textAnchor="middle">{v}</text>
          </g>;
        })}
        <text x={cx} y={bottom + 26} fontSize="13" fill={SOFT} textAnchor="middle">{item.label}</text>
      </g>;
    })}
  </svg>;
}

export function Waffles({ panels, svgRef }: {
  panels: Array<{ total: number; highlighted: number; title: string }>; svgRef?: RefObject<SVGSVGElement | null>;
}) {
  const perRow = 40, r = 3.9, step = 10.8, blockW = perRow * step, gapX = 64;
  const rowsMax = Math.max(...panels.map((p) => Math.ceil(p.total / perRow)));
  const width = panels.length * blockW + (panels.length - 1) * gapX + 20, height = 70 + rowsMax * step + 20;
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg"
    aria-label={panels.map((p) => `${p.highlighted} of ${p.total} ${p.title}`).join("; ")} style={{ maxWidth: width }}>
    {panels.map((p, k) => {
      const x0 = 10 + k * (blockW + gapX);
      return <g key={p.title}>
        <text x={x0} y="22" fontSize="15" fontWeight="600" fill={INK}>{p.title}</text>
        <text x={x0} y="44" fontSize="13" fill={SOFT}>{p.highlighted} of {p.total} metabolites significant</text>
        {Array.from({ length: p.total }, (_, i) => <circle key={i} cx={x0 + (i % perRow) * step + r} cy={66 + Math.floor(i / perRow) * step}
          r={r} fill={i < p.highlighted ? COLOURS.warm : "#e6e1da"} />)}
      </g>;
    })}
  </svg>;
}

export function Scatter({ points, xLabel, yLabel, svgRef }: {
  points: Array<{ x: number; y: number; label: string }>; xLabel: string; yLabel: string; svgRef?: RefObject<SVGSVGElement | null>;
}) {
  const width = 900, height = 420, left = 80, right = 30, top = 30, bottom = 360;
  const max = Math.max(0.1, ...points.flatMap((p) => [p.x, p.y])) * 1.08;
  const sx = (v: number) => left + (v / max) * (width - left - right);
  const sy = (v: number) => bottom - (v / max) * (bottom - top);
  const ticks = [0, max / 4, max / 2, (3 * max) / 4].map((t) => Number(t.toFixed(1)));
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg" aria-label={`${yLabel} against ${xLabel}`} style={{ maxWidth: 820 }}>
    {ticks.map((t) => <g key={t}>
      <line x1={left} x2={width - right} y1={sy(t)} y2={sy(t)} stroke="#efebe5" />
      <text x={left - 10} y={sy(t) + 4} fontSize="12" fill={FAINT} textAnchor="end">{t}</text>
      <text x={sx(t)} y={bottom + 20} fontSize="12" fill={FAINT} textAnchor="middle">{t}</text>
    </g>)}
    <line x1={sx(0)} y1={sy(0)} x2={sx(max)} y2={sy(max)} stroke="#bdb6ae" strokeDasharray="5 5" />
    <text x={sx(max * 0.72)} y={sy(max * 0.72) + 22} fontSize="12.5" fill={FAINT} transform={`rotate(-${Math.atan((bottom - top) / (width - left - right)) * 180 / Math.PI}, ${sx(max * 0.72)}, ${sy(max * 0.72) + 22})`}>equal to healthy drift</text>
    {points.map((p) => <circle key={p.label} cx={sx(p.x)} cy={sy(p.y)} r="7" fill={COLOURS.warm} stroke="#fff" strokeWidth="2">
      <title>{`${p.label}: PAH difference ${p.y.toFixed(2)}, largest healthy drift ${p.x.toFixed(2)} (log2)`}</title></circle>)}
    <text x={(left + width - right) / 2} y={height - 14} fontSize="13" fill={SOFT} textAnchor="middle">{xLabel}</text>
    <text x="18" y={(top + bottom) / 2} fontSize="13" fill={SOFT} textAnchor="middle" transform={`rotate(-90, 18, ${(top + bottom) / 2})`}>{yLabel}</text>
  </svg>;
}

/* Dot plot (pathway × time): colour = signed score, dot area grows as q shrinks, dark ring = significant. */
export type DotCell = { value: number | null; q: number | null; pass: boolean | null; title: string } | null;

export function DotPlot({ rows, groups, cell, labels, svgRef, rowLabelWidth = 220 }: {
  rows: string[]; groups: Array<{ title?: string; cols: Col[] }>; cell: (row: string, col: string) => DotCell;
  labels: [string, string]; svgRef?: RefObject<SVGSVGElement | null>; rowLabelWidth?: number;
}) {
  const nCols = groups.reduce((n, g) => n + g.cols.length, 0);
  const groupGap = 26, colW = Math.min(76, (940 - rowLabelWidth - groupGap * (groups.length - 1)) / Math.max(nCols, 1)), rowH = 44;
  const top = 64, width = rowLabelWidth + nCols * colW + groupGap * (groups.length - 1) + 12, height = top + rows.length * rowH + 78;
  const values = rows.flatMap((r) => groups.flatMap((g) => g.cols.map((c) => Math.abs(cell(r, c.key)?.value ?? 0))));
  const scale = Math.max(1, ...values);
  const radius = (q: number | null) => 4 + Math.min(1, Math.max(0, -Math.log10(Math.max(q ?? 1, 1e-12))) / 6) * 13;
  let x0 = rowLabelWidth;
  const placed = groups.map((g) => { const start = x0; x0 += g.cols.length * colW + groupGap; return { ...g, start }; });
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg"
    aria-label={`Dot plot, ${rows.length} rows by ${nCols} columns`} style={{ maxWidth: width }}>
    {placed.map((g) => <g key={g.start}>
      {g.title && <text x={g.start + 4} y="18" fontSize="14" fontWeight="600" fill={INK}>{g.title}</text>}
      {g.cols.map((c, i) => <text key={c.key} x={g.start + i * colW + colW / 2} y={top - 16} fontSize="12" fill={SOFT} textAnchor="middle">{c.label}</text>)}
      {rows.map((r, j) => g.cols.map((c, i) => {
        const v = cell(r, c.key), cx = g.start + i * colW + colW / 2, cy = top + j * rowH + rowH / 2;
        return <g key={r + c.key}>
          <title>{v?.title ?? `${r} · ${c.label}: not measured`}</title>
          <line x1={cx - colW / 2 + 4} x2={cx + colW / 2 - 4} y1={cy} y2={cy} stroke="#f1ede7" />
          {v ? <circle cx={cx} cy={cy} r={radius(v.q)} fill={diverging(v.value, scale)} stroke={v.pass ? INK : "#fff"} strokeWidth={v.pass ? 2 : 1.5} />
            : <circle cx={cx} cy={cy} r="3" fill="none" stroke={RULE} />}
        </g>;
      }))}
    </g>)}
    {rows.map((r, j) => <text key={r} x={rowLabelWidth - 14} y={top + j * rowH + rowH / 2 + 5} fontSize="14" fill={INK} textAnchor="end">{r}</text>)}
    <Legend x={rowLabelWidth + 110} y={height - 50} scale={scale} labels={labels} dotKey={false} />
    <g transform={`translate(${rowLabelWidth + 110}, ${height - 18})`}>
      <text x="0" y="4" fontSize="13" fill={SOFT} textAnchor="end">q</text>
      {[0.5, 0.05, 0.001].map((q, i) => <g key={q}><circle cx={18 + i * 70} cy="0" r={radius(q)} fill="#e6e1da" /><text x={18 + i * 70 + radius(q) + 5} y="4" fontSize="12" fill={SOFT}>{q}</text></g>)}
      <circle cx={250} cy="0" r="7" fill="#fff" stroke={INK} strokeWidth="2" /><text x={262} y="4" fontSize="12.5" fill={SOFT}>ring: significant</text>
    </g>
  </svg>;
}

/* Small multiples: one panel per molecule, value over time, one line per series. */
export type Series = { name: string; colour: string; points: Array<{ x: number; y: number | null; pass: boolean | null; title: string }> };

export function SmallMultiples({ panels, xLabels, yLabel, svgRef, perRow = 3 }: {
  panels: Array<{ title: string; series: Series[] }>; xLabels: string[]; yLabel: string;
  svgRef?: RefObject<SVGSVGElement | null>; perRow?: number;
}) {
  const pw = 290, ph = 170, gapX = 22, gapY = 34, top = 44, padL = 40, padB = 26;
  const rows = Math.ceil(panels.length / perRow);
  const width = perRow * pw + (perRow - 1) * gapX + 10, height = top + rows * (ph + gapY) + 10;
  const ys = panels.flatMap((p) => p.series.flatMap((s) => s.points.map((pt) => pt.y ?? 0)));
  const lim = Math.max(0.1, ...ys.map(Math.abs)) * 1.1;
  const names = Array.from(new Map(panels.flatMap((p) => p.series.map((s) => [s.name, s.colour] as const))).entries());
  return <svg ref={svgRef} viewBox={`0 0 ${width} ${height}`} width="100%" role="img" className="lr-svg" aria-label={`${yLabel} over time, ${panels.length} panels`} style={{ maxWidth: width }}>
    {names.map(([name, colour], i) => <g key={name}><line x1={10 + i * 150} x2={32 + i * 150} y1="16" y2="16" stroke={colour} strokeWidth="2.5" />
      <circle cx={21 + i * 150} cy="16" r="4" fill={colour} /><text x={40 + i * 150} y="20" fontSize="13" fill={INK}>{name}</text></g>)}
    <text x={width - 10} y="20" fontSize="12" fill={FAINT} textAnchor="end">{yLabel} · filled point: significant · shared y axis ±{lim.toFixed(2)}</text>
    {panels.map((p, k) => {
      const x0 = (k % perRow) * (pw + gapX), y0 = top + Math.floor(k / perRow) * (ph + gapY);
      const sx = (i: number) => x0 + padL + (xLabels.length > 1 ? i / (xLabels.length - 1) : 0.5) * (pw - padL - 12);
      const sy = (v: number) => y0 + 18 + (1 - (v + lim) / (2 * lim)) * (ph - 18 - padB);
      return <g key={p.title}>
        <rect x={x0} y={y0} width={pw} height={ph} rx="12" fill="#faf8f5" />
        <text x={x0 + 12} y={y0 + 14} fontSize="13" fontWeight="600" fill={INK}>{p.title}</text>
        <line x1={x0 + padL} x2={x0 + pw - 12} y1={sy(0)} y2={sy(0)} stroke="#bdb6ae" strokeDasharray="3 3" />
        {xLabels.map((l, i) => <text key={l} x={sx(i)} y={y0 + ph - 8} fontSize="11" fill={SOFT} textAnchor="middle">{l}</text>)}
        {p.series.map((s) => {
          const pts = s.points.filter((pt) => pt.y !== null);
          return <g key={s.name}>
            <polyline points={pts.map((pt) => `${sx(pt.x)},${sy(pt.y!)}`).join(" ")} fill="none" stroke={s.colour} strokeWidth="2" strokeLinejoin="round" />
            {pts.map((pt) => <circle key={pt.x} cx={sx(pt.x)} cy={sy(pt.y!)} r="4.5" fill={pt.pass ? s.colour : "#fff"} stroke={s.colour} strokeWidth="2"><title>{pt.title}</title></circle>)}
          </g>;
        })}
      </g>;
    })}
  </svg>;
}
