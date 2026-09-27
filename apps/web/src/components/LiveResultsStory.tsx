import { useEffect, useMemo, useRef, useState, type ReactNode, type RefObject } from "react";
import {
  api, ApiError, type AnalysisResponse, type ColumnResult, type MappingPreview, type MetabolitePreview,
  type StoryPanel, type StoryRecord, type StoryResponse,
} from "../api/client";
import {
  ALL, AGREEMENT, COMBINED_SEX, DEFAULT_FILTERS, DIMENSIONS, agreementClass, columnContext, columnDirection,
  columnSignificance, filterRecords, matchesFilters, optionAvailable, optionLabel, optionValues,
  significancePass, significanceValue, unavailableReason, unavailableRuleReason,
  type Dimension, type FilterContext, type LiveFilters, type Significance,
} from "./liveDashboardModel";
import "./LiveDashboard.css";

// ---- constants ----------------------------------------------------------------------------------

type LoadState<T> = { kind: "idle" } | { kind: "loading" } | { kind: "error"; message: string } | { kind: "ready"; data: T };
type StepId = "1" | "2" | "3" | "4" | "5" | "data" | "catalog";

const STEPS: Array<{ id: StepId; title: string; subtitle: string; panels: string[] }> = [
  { id: "1", title: "Human protein", subtitle: "Nine PAH markers", panels: ["human_muscle_protein"] },
  { id: "2", title: "RNA across tissues", subtitle: "The discordance event", panels: [
    "human_muscle_rna", "muscle_oxphos", "blood_ranks", "blood_go", "blood_biocarta"] },
  { id: "3", title: "Rat training", subtitle: "Weeks, not hours", panels: ["rat_muscle_protein"] },
  { id: "4", title: "Metabolites", subtitle: "Patients and sampling", panels: [
    "metabolomics_exercise", "metabolomics_rest", "metabolomics_context"] },
  { id: "5", title: "Second disease", subtitle: "Type 2 diabetes", panels: [] },
  { id: "data", title: "Your data", subtitle: "Genes, proteins, metabolites", panels: [] },
  { id: "catalog", title: "Evidence catalog", subtitle: "Molecule × layer × time", panels: [] },
];

const FILTER_LABELS: Record<Dimension, string> = {
  exercise: "Exercise", contrast: "Contrast", timepoint: "Time", layer: "Layer",
  species: "Species", sex: "Sex", tissue: "Tissue",
};
const RULE_LABEL: Record<Significance, string> = { raw: "Raw P < 0.05", bh: "BH q < 0.05", bonferroni: "Bonferroni < 0.05" };

const PANEL_CONTEXT: Record<string, { rest: string; exercise: string }> = {
  human_muscle_protein: {
    rest: "Malenfant 2015: nine proteins lower in resting PAH vastus lateralis than in controls (4 vs 4).",
    exercise: "MoTrPAC: sedentary healthy adults after one endurance bout; muscle protein, EE-CON.",
  },
  human_muscle_rna: {
    rest: "Disease direction from Malenfant's protein result; no PAH RNA enters this panel.",
    exercise: "MoTrPAC: healthy muscle RNA after one endurance bout, EE-CON.",
  },
  muscle_oxphos: {
    rest: "No PAH test enters this panel.",
    exercise: "MoTrPAC: healthy muscle RNA, GO oxidative phosphorylation (CAMERA).",
  },
  blood_ranks: {
    rest: "GSE33463: all genes ranked by IPAH vs healthy PBMC t (30 vs 41).",
    exercise: "MoTrPAC: healthy whole-blood RNA ranked by endurance vs control z.",
  },
  blood_go: {
    rest: "GSE33463: IPAH vs healthy PBMC GO pathway tests.",
    exercise: "MoTrPAC: healthy whole-blood GO tests; labels joined by exact name.",
  },
  blood_biocarta: {
    rest: "Cheadle 2012: six BioCarta pathways lower in IPAH blood cells (PAGE z).",
    exercise: "MoTrPAC: healthy whole-blood CAMERA z; modality shown per row.",
  },
  rat_muscle_protein: {
    rest: "Disease direction from Malenfant's human PAH muscle protein result.",
    exercise: "MoTrPAC: healthy rats trained 1, 2, 4 or 8 weeks vs sedentary rats.",
  },
  metabolomics_exercise: {
    rest: "ST000763: plasma at rest in SSc-PAH, SSc without PH and healthy people.",
    exercise: "ST000763: the same people at peak exercise; rest/peak pairs are inferred.",
  },
  metabolomics_rest: {
    rest: "ST000763: resting PAH vs healthy, and vs normal-pressure SSc sampled in the same catheter setting.",
    exercise: "No exercise data in this test; MoTrPAC context follows below.",
  },
  metabolomics_context: {
    rest: "ST000763: the 41 resting PAH-vs-healthy metabolites.",
    exercise: "MoTrPAC: healthy blood — exercise effect (EE-CON), exercisers only (EE-EE), control drift (CON-CON).",
  },
};

const T2D_EXAMPLE = "type2_diabetes_muscle_mootha2003";
const PAH_EXAMPLE = "pah_muscle_lower9_malenfant2015";
const GENE_TEMPLATE = "gene_symbol,direction\nNDUFA9,-1\nUQCRC1,-1\nATP5F1B,-1\n";
const METABOLITE_COLUMNS = ["list_id", "source", "source_type", "comparison", "metabolite_name", "refmet_name",
  "pah_log2_effect", "pah_direction", "pah_q_value", "is_reported_hit", "caveat"];
const METABOLITE_TEMPLATE = METABOLITE_COLUMNS.join(",") +
  "\nmy_list,Author 2026,paper_reported,PAH vs control at rest,Succinate,Succinic acid,0.6,up,0.01,True,\n";

// ---- small helpers ------------------------------------------------------------------------------

function fmt(value: number | null | undefined): string {
  if (value === null || value === undefined || !Number.isFinite(value)) return "—";
  if (value === 0) return "0";
  if (Math.abs(value) < 0.001) return value.toExponential(2);
  return value.toFixed(3).replace(/0+$/, "").replace(/\.$/, "");
}
function errorText(error: unknown): string {
  return error instanceof ApiError ? `${error.status}: ${error.detail}` : String(error);
}
function readFilters(): LiveFilters {
  const params = new URLSearchParams(window.location.search);
  const next = { ...DEFAULT_FILTERS };
  for (const key of DIMENSIONS) next[key] = params.get("lr_" + key) ?? DEFAULT_FILTERS[key];
  const rule = params.get("lr_significance");
  if (rule === "raw" || rule === "bh" || rule === "bonferroni") next.significance = rule;
  return next;
}
function readStep(): StepId {
  const selected = new URLSearchParams(window.location.search).get("lr_step");
  return STEPS.find((step) => step.id === selected)?.id ?? "1";
}
function writeUrl(filters: LiveFilters, step: StepId): void {
  const url = new URL(window.location.href);
  for (const key of DIMENSIONS) {
    if (filters[key] === ALL) url.searchParams.delete("lr_" + key);
    else url.searchParams.set("lr_" + key, filters[key]);
  }
  url.searchParams.set("lr_significance", filters.significance);
  url.searchParams.set("lr_step", step);
  window.history.replaceState(null, "", url);
}
function download(name: string, contents: Blob): void {
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
function csvDownload(name: string, rows: Array<Record<string, unknown>>): void {
  const columns = rows.length ? Object.keys(rows[0]) : ["empty"];
  const contents = [columns.map(csvCell).join(","), ...rows.map((row) =>
    columns.map((column) => csvCell(row[column])).join(","))].join("\r\n");
  download(name, new Blob([contents], { type: "text/csv;charset=utf-8" }));
}
async function svgPng(svg: SVGSVGElement | null, name: string): Promise<void> {
  if (!svg) return;
  const copy = svg.cloneNode(true) as SVGSVGElement;
  copy.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)], { type: "image/svg+xml;charset=utf-8" }));
  const image = new Image();
  try {
    await new Promise<void>((resolve, reject) => {
      image.onload = () => resolve();
      image.onerror = () => reject(new Error("Could not render the figure."));
      image.src = url;
    });
    const scale = 2;
    const canvas = document.createElement("canvas");
    canvas.width = (Number(svg.getAttribute("width")) || 900) * scale;
    canvas.height = (Number(svg.getAttribute("height")) || 400) * scale;
    const context = canvas.getContext("2d");
    if (!context) throw new Error("Canvas is unavailable.");
    context.scale(scale, scale);
    context.fillStyle = "#fbfaf8";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0);
    const blob = await new Promise<Blob | null>((resolve) => canvas.toBlob(resolve, "image/png"));
    if (!blob) throw new Error("PNG export failed.");
    download(name, blob);
  } finally { URL.revokeObjectURL(url); }
}
function readFile(file: File | undefined, onText: (text: string) => void): void {
  if (!file) return;
  file.text().then(onText).catch(() => undefined);
}

function recordRow(record: StoryRecord, rule: Significance): Record<string, unknown> {
  return {
    panel: record.panel_id, record_type: record.record_type, molecule: record.molecule ?? "", label: record.label,
    species: record.species, tissue: record.tissue, layer: record.layer, exercise: record.exercise,
    contrast: record.contrast, timepoint: record.timepoint, sex: record.sex,
    disease_contrast: record.disease_contrast, disease_direction: record.disease_direction,
    disease_value: record.disease_value, disease_unit: record.disease_unit, disease_q: record.disease_q_value,
    motrpac_measured: record.motrpac_measured, motrpac_direction: record.motrpac_direction,
    estimate: record.value, unit: record.unit, n: record.n, raw_p: record.p_value, bh_q: record.q_value,
    bonferroni: record.bonferroni_p, bh_family: record.correction_family,
    rule, passes_rule: significancePass(record, rule), agreement: agreementClass(record, rule),
    source: record.source,
  };
}

// ---- filter bar -------------------------------------------------------------------------------

function FilterBar({ contexts, filters, onChange, onReset }: {
  contexts: FilterContext[]; filters: LiveFilters; onChange: (next: LiveFilters) => void; onReset: () => void;
}) {
  const exact = contexts.filter((context) => matchesFilters(context, filters)).length;
  const active = DIMENSIONS.filter((dimension) => filters[dimension] !== ALL).length;
  return <div className="lr-filters" aria-label="Filters for every step">
    <div className="lr-filter-heading">
      <div><strong>Filters</strong><span>{exact} result contexts match · applied to every step</span></div>
      {active > 0 && <button type="button" className="lr-text-button" onClick={onReset}>Reset</button>}
    </div>
    <div className="lr-filter-grid">
      {DIMENSIONS.map((dimension) => <details className="lr-filter" key={dimension}>
        <summary><span>{FILTER_LABELS[dimension]}</span><strong>{filters[dimension] === ALL ? "All" : optionLabel(filters[dimension])}</strong></summary>
        <div className="lr-filter-menu">{optionValues(contexts, dimension).map((option) => {
          const available = optionAvailable(contexts, filters, dimension, option);
          const selected = filters[dimension] === option;
          return <button type="button" key={option}
            title={available ? undefined : unavailableReason(dimension, option, contexts, filters)}
            aria-disabled={!available && !selected} aria-pressed={selected}
            className={"lr-option" + (!available ? " lr-option--unavailable" : "") + (selected ? " lr-option--selected" : "")}
            onClick={(event) => {
              if (!available && !selected) return;
              onChange({ ...filters, [dimension]: option });
              event.currentTarget.closest("details")?.removeAttribute("open");
            }}>
            {option === ALL ? "All" : optionLabel(option)}
            {!available && <small>{unavailableReason(dimension, option, contexts, filters)}</small>}
          </button>;
        })}</div>
      </details>)}
      <details className="lr-filter">
        <summary><span>Significance</span><strong>{RULE_LABEL[filters.significance]}</strong></summary>
        <div className="lr-filter-menu">{(Object.keys(RULE_LABEL) as Significance[]).map((value) =>
          <button type="button" key={value} aria-pressed={filters.significance === value}
            className={"lr-option" + (filters.significance === value ? " lr-option--selected" : "")}
            onClick={(event) => { onChange({ ...filters, significance: value }); event.currentTarget.closest("details")?.removeAttribute("open"); }}>
            {RULE_LABEL[value]}
            <small>{value === "raw" ? "Unadjusted" : value === "bh" ? "The source study's own BH family" : "Same family as BH; only where it can be verified"}</small>
          </button>)}
        </div>
      </details>
    </div>
  </div>;
}

// ---- figure and table -------------------------------------------------------------------------

function plottedRecords(panel: StoryPanel, records: StoryRecord[], filters: LiveFilters): StoryRecord[] {
  let selected = records.filter((record) => record.value !== null);
  if (panel.id === "blood_go") selected = selected.filter((record) => record.record_type === "pathway_time_summary");
  if (panel.id === "rat_muscle_protein") {
    const official = selected.filter((record) => record.record_type === "official_sex_specific");
    if (official.length) selected = official;
  }
  if (panel.id === "metabolomics_context" && filters.contrast === ALL && filters.timepoint === ALL) {
    selected = selected.filter((record) => record.record_type === "metabolite_context");
  }
  return selected.slice(0, 40);
}

function Plot({ records, rule, svgRef }: { records: StoryRecord[]; rule: Significance; svgRef: RefObject<SVGSVGElement | null> }) {
  const width = 900, rowHeight = 24, top = 84, origin = 470, span = 290;
  const height = Math.max(150, top + records.length * rowHeight + 22);
  const scale = Math.max(1e-9, ...records.map((record) => Math.abs(record.value ?? 0)));
  const units = Array.from(new Set(records.map((record) => record.unit).filter(Boolean)));
  return <svg ref={svgRef} className="lr-plot" width={width} height={height} viewBox={`0 0 ${width} ${height}`}
    role="img" aria-label="Estimates for the displayed rows">
    <rect width={width} height={height} fill="#fbfaf8" />
    <text x="16" y="24" fontSize="13" fontWeight="600" fill="#33312e">{units.join(" · ") || "Estimate"}</text>
    <text x="16" y="42" fontSize="11" fill="#6d6862">Filled: passes {RULE_LABEL[rule]} · hollow: rule not checkable · grey: does not pass</text>
    <line x1={origin} x2={origin} y1={top - 12} y2={height - 12} stroke="#bdb6ae" />
    <text x={origin - span} y={top - 20} fontSize="10" fill="#9a938b">− {fmt(scale)}</text>
    <text x={origin + span} y={top - 20} fontSize="10" fill="#9a938b" textAnchor="end">+ {fmt(scale)}</text>
    <text x={origin} y={top - 20} fontSize="10" fill="#9a938b" textAnchor="middle">0</text>
    {records.map((record, index) => {
      const y = top + index * rowHeight;
      const x = origin + ((record.value ?? 0) / scale) * span;
      const pass = significancePass(record, rule);
      const label = record.molecule
        ? record.molecule + (record.timepoint ? " · " + optionLabel(record.timepoint) : "") +
          (record.sex && !["all", "both"].includes(record.sex) ? " · " + record.sex : "") +
          (record.record_type === "metabolite_timecourse" ? " · " + record.contrast : "")
        : record.timepoint ? optionLabel(record.timepoint) + (record.contrast?.includes("CON") ? "" : " · " + record.label) : record.label;
      const short = label.length > 58 ? label.slice(0, 55) + "…" : label;
      const colour = pass ? "#8c2740" : "#bdb6ae";
      return <g key={record.id}>
        <title>{`${label}\n${fmt(record.value)} ${record.unit ?? ""}\nraw P ${fmt(record.p_value)} · BH q ${fmt(record.q_value)} · Bonferroni ${fmt(record.bonferroni_p)}`}</title>
        <text x="16" y={y + 4} fontSize="11" fill="#33312e">{short}</text>
        <line x1={origin} x2={x} y1={y} y2={y} stroke={colour} strokeWidth="2" />
        <circle cx={x} cy={y} r={5.5} fill={pass === null ? "#fffefd" : colour} stroke={pass === null ? "#8c2740" : "#fbfaf8"} strokeWidth="2" />
        <text x={width - 16} y={y + 4} fontSize="11" fill="#6d6862" textAnchor="end">{fmt(record.value)}</text>
      </g>;
    })}
    {!records.length && <text x="16" y={top + 20} fontSize="12" fill="#6d6862">No numeric estimate for this selection.</text>}
  </svg>;
}

function RuleCell({ record, rule }: { record: StoryRecord; rule: Significance }) {
  const pass = significancePass(record, rule);
  if (pass === null) return <td className="lr-muted" title={unavailableRuleReason(record, rule)}>not checkable</td>;
  return <td className={pass ? "lr-pass" : "lr-muted"}>{pass ? "passes" : "does not pass"}</td>;
}

function RecordTable({ records, rule }: { records: StoryRecord[]; rule: Significance }) {
  return <div className="lr-table-wrap"><table className="lr-table">
    <thead><tr><th>Result</th><th>Context</th><th>Disease direction</th><th>MoTrPAC direction</th>
      <th>Estimate</th><th>Raw P</th><th>BH q</th><th>Bonferroni</th><th>{RULE_LABEL[rule]}</th></tr></thead>
    <tbody>{records.map((record) => <tr key={record.id}>
      <th scope="row">{record.molecule ?? record.label}<small>{record.molecule ? record.label : record.record_type.replace(/_/g, " ")}</small></th>
      <td>{record.species} · {record.tissue} · {record.layer}<small>{optionLabel(record.contrast ?? "")} · {optionLabel(record.timepoint ?? "")}{record.sex && !["all", "both"].includes(record.sex) ? " · " + record.sex : ""}</small></td>
      <td>{record.disease_direction ?? "—"}{record.disease_value !== null && record.disease_value !== undefined && <small>{fmt(record.disease_value)} {record.disease_unit ?? ""}</small>}</td>
      <td>{record.motrpac_measured === false ? "not measured" : record.motrpac_direction ?? "—"}</td>
      <td>{fmt(record.value)}<small>{record.unit ?? ""}{record.n ? ` · n ${record.n}` : ""}</small></td>
      <td>{fmt(record.p_value)}</td><td>{fmt(record.q_value)}</td><td>{fmt(record.bonferroni_p)}</td>
      <RuleCell record={record} rule={rule} />
    </tr>)}</tbody>
  </table></div>;
}

function FigureTools({ view, setView, onCsv, onPng, note }: {
  view: "plot" | "table"; setView: (view: "plot" | "table") => void; onCsv: () => void; onPng?: () => void; note?: ReactNode;
}) {
  return <div className="lr-figure-tools">
    <div className="lr-segment" role="group" aria-label="View">
      <button type="button" aria-pressed={view === "plot"} onClick={() => setView("plot")}>Figure</button>
      <button type="button" aria-pressed={view === "table"} onClick={() => setView("table")}>Table</button>
    </div>
    <span>{note}</span>
    <button type="button" className="lr-outline-button" onClick={onCsv}>CSV</button>
    {onPng && <button type="button" className="lr-outline-button" onClick={onPng}>PNG</button>}
  </div>;
}

// ---- one story panel --------------------------------------------------------------------------

function PanelCard({ panel, filters, onUseContext }: { panel: StoryPanel; filters: LiveFilters; onUseContext: (context: FilterContext) => void }) {
  const [view, setView] = useState<"plot" | "table">("plot");
  const [exportError, setExportError] = useState<string | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);
  const records = useMemo(() => filterRecords(panel.records, filters), [panel.records, filters]);
  const plotted = useMemo(() => plottedRecords(panel, records, filters), [panel, records, filters]);
  const context = PANEL_CONTEXT[panel.id];
  const passing = records.filter((record) => significancePass(record, filters.significance) === true).length;
  const checkable = records.filter((record) => significancePass(record, filters.significance) !== null).length;
  const families = Array.from(new Set(records.map((record) => record.correction_family).filter(Boolean)));
  return <article className="lr-panel" id={"lr-" + panel.id}>
    <header className="lr-panel-head">
      <span className="lr-eyebrow">{panel.title}</span>
      <h3>{panel.question}</h3>
    </header>
    <div className="lr-contrast-pair">
      <div><span>Rest · disease study</span><p>{context?.rest}</p></div>
      <div><span>Exercise · healthy MoTrPAC</span><p>{context?.exercise}</p></div>
    </div>
    <div className="lr-result"><span>Result</span><strong>{panel.observation}</strong></div>
    <p className="lr-caveat"><span>Caveat</span>{panel.caveat}</p>
    <details className="lr-interpretation"><summary>Interpretation and hypothesis</summary>
      <p><strong>Interpretation:</strong> {panel.interpretation}</p>
      <p><strong>Hypothesis:</strong> {panel.hypothesis}</p>
    </details>
    {records.length ? <div className="lr-figure">
      <FigureTools view={view} setView={setView}
        note={checkable ? `${passing} of ${checkable} checkable rows pass ${RULE_LABEL[filters.significance]}` : `${RULE_LABEL[filters.significance]} not checkable here`}
        onCsv={() => csvDownload(panel.id + ".csv", records.map((record) => recordRow(record, filters.significance)))}
        onPng={() => svgPng(svgRef.current, panel.id + ".png").catch((cause) => setExportError(errorText(cause)))} />
      <div className={view === "plot" ? "lr-chart-scroll" : "lr-hidden"}>
        <Plot records={plotted} rule={filters.significance} svgRef={svgRef} />
      </div>
      {view === "table" && <RecordTable records={records} rule={filters.significance} />}
      {view === "plot" && plotted.length < records.filter((record) => record.value !== null).length &&
        <p className="lr-small">The figure shows {plotted.length} rows; the table and CSV hold all {records.length}. Narrow the filters to plot others.</p>}
      {exportError && <p className="lr-error">{exportError}</p>}
      <details className="lr-sources"><summary>Sources and correction families</summary>
        {Array.from(new Set(records.map((record) => record.source))).map((source) => <p key={source}><code>MoTrPAC Hackathon/{source}</code></p>)}
        {families.map((family) => <p key={family}>BH family: {family}</p>)}
      </details>
    </div> : <div className="lr-empty">
      <strong>No row in this panel matches the filters.</strong>
      <p>The result above still describes the whole panel.</p>
      <button type="button" className="lr-outline-button" onClick={() => onUseContext(panel.records[0])}>Show this panel's data</button>
    </div>}
  </article>;
}

// ---- step 2: RNA versus protein in the same muscle ---------------------------------------------

type LayerPair = { gene: string; timepoint: string; rna?: StoryRecord; protein?: StoryRecord; label: string };

function pairClass(rna: StoryRecord | undefined, protein: StoryRecord | undefined, rule: Significance): string {
  if (!rna || !protein) return "one layer not measured";
  const r = significancePass(rna, rule), p = significancePass(protein, rule);
  if (r === null || p === null) return "rule not checkable";
  if (r && p) return (rna.value ?? 0) * (protein.value ?? 0) > 0 ? "concordant change" : "discordant change";
  if (r) return "discordant: RNA only";
  if (p) return "discordant: protein only";
  return "concordant: neither changes";
}

function layerPairs(story: StoryResponse, rule: Significance, filters: LiveFilters): LayerPair[] {
  const panel = (id: string) => story.panels.find((item) => item.id === id)?.records ?? [];
  const key = (record: StoryRecord) => `${String(record.metrics?.uniprot)}|${record.timepoint}`;
  const proteins = new Map(panel("human_muscle_protein").map((record) => [key(record), record]));
  return panel("human_muscle_rna")
    .filter((record) => filters.timepoint === ALL || record.timepoint === filters.timepoint)
    .map((rna) => {
      const protein = proteins.get(key(rna));
      return { gene: rna.molecule ?? "", timepoint: rna.timepoint ?? "", rna, protein, label: pairClass(rna, protein, rule) };
    });
}

function humanMuscleVisible(filters: LiveFilters): boolean {
  const ok = (dimension: Dimension, values: string[]) => filters[dimension] === ALL || values.includes(filters[dimension]);
  return ok("species", ["human"]) && ok("tissue", ["skeletal muscle"]) && ok("contrast", ["EE-CON"]) &&
    ok("exercise", ["endurance"]) && ok("sex", [COMBINED_SEX]);
}

function LayerDiscordance({ story, filters }: { story: StoryResponse; filters: LiveFilters }) {
  const pairs = layerPairs(story, filters.significance, filters);
  const times = Array.from(new Set(pairs.map((pair) => pair.timepoint)));
  const genes = Array.from(new Set(pairs.map((pair) => pair.gene)));
  if (!humanMuscleVisible(filters) || !pairs.length) return null;
  const tone = (label: string) => label.startsWith("concordant change") ? "same" : label.startsWith("discordant") ? "discordant" : "quiet";
  return <article className="lr-panel lr-panel--feature" id="lr-layer-discordance">
    <header className="lr-panel-head">
      <span className="lr-eyebrow">The discordance event</span>
      <h3>Do RNA and protein for the nine genes change together in the same healthy muscle?</h3>
    </header>
    <div className="lr-result"><span>Result</span><strong>
      {times.map((time) => {
        const at = pairs.filter((pair) => pair.timepoint === time);
        const rnaOnly = at.filter((pair) => pair.label === "discordant: RNA only").length;
        return `${optionLabel(time)}: ${rnaOnly}/${at.length} RNA-only`;
      }).join(" · ")}
    </strong></div>
    <p className="lr-caveat"><span>Caveat</span>Same MoTrPAC muscle biopsies, separate assays. A protein change could follow later than 24 h, the last sampled time.</p>
    <div className="lr-matrix" style={{ gridTemplateColumns: `110px repeat(${times.length}, minmax(150px, 1fr))` }}>
      <span />
      {times.map((time) => <strong key={time}>{optionLabel(time)}</strong>)}
      {genes.map((gene) => <FragmentRow key={gene} gene={gene} times={times} pairs={pairs} tone={tone} rule={filters.significance} />)}
    </div>
    <p className="lr-small">Each cell: RNA and protein EE-CON log2 fold change, marked • where {RULE_LABEL[filters.significance]} passes.</p>
    <div className="lr-figure-tools"><span />
      <button type="button" className="lr-outline-button" onClick={() => csvDownload("rna_protein_pairs.csv", pairs.map((pair) => ({
        gene: pair.gene, timepoint: pair.timepoint, rna_logFC: pair.rna?.value, rna_q: pair.rna?.q_value,
        rna_bonferroni: pair.rna?.bonferroni_p, protein_logFC: pair.protein?.value, protein_q: pair.protein?.q_value,
        protein_bonferroni: pair.protein?.bonferroni_p, rule: filters.significance, class: pair.label,
      })))}>CSV</button></div>
  </article>;
}

function FragmentRow({ gene, times, pairs, tone, rule }: {
  gene: string; times: string[]; pairs: LayerPair[]; tone: (label: string) => string; rule: Significance;
}) {
  const mark = (record: StoryRecord | undefined) => record && significancePass(record, rule) ? " •" : "";
  return <>
    <span className="lr-matrix-gene">{gene}</span>
    {times.map((time) => {
      const pair = pairs.find((item) => item.gene === gene && item.timepoint === time);
      if (!pair) return <span key={time} className="lr-cell lr-cell--quiet">—</span>;
      return <span key={time} className={"lr-cell lr-cell--" + tone(pair.label)} title={pair.label}>
        <small>RNA {fmt(pair.rna?.value)}{mark(pair.rna)}</small>
        <small>Protein {fmt(pair.protein?.value)}{mark(pair.protein)}</small>
        <em>{pair.label}</em>
      </span>;
    })}
  </>;
}

// ---- engine comparison tables (step 5 and "your data") --------------------------------------

function ColumnTable({ runs, filters }: { runs: Array<{ name: string; data: AnalysisResponse }>; filters: LiveFilters }) {
  const [view, setView] = useState<"plot" | "table">("table");
  const base = runs[0].data.columns.filter((column) => matchesFilters(columnContext(column), filters));
  const byId = runs.map((run) => new Map(run.data.columns.map((column) => [column.column_id, column])));
  const counts = runs.map((run, index) => {
    const shown = base.map((column) => byId[index].get(column.column_id)).filter(Boolean) as ColumnResult[];
    const tally = (label: string) => shown.filter((column) => columnDirection(column, filters.significance) === label).length;
    return { name: run.name, opposite: tally("opposite direction to disease"), same: tally("same direction as disease"), total: shown.length };
  });
  const rows = base.map((column) => ({
    comparison: column.comparison_label,
    ...Object.fromEntries(runs.flatMap((run, index) => {
      const hit = byId[index].get(column.column_id);
      return [
        [`${run.name}_measured`, hit ? `${hit.n_measured}` : ""],
        [`${run.name}_opposite_members`, hit?.n_opposed ?? ""],
        [`${run.name}_camera_t`, hit?.camera_t ?? ""],
        [`${run.name}_${filters.significance}`, hit ? columnSignificance(hit, filters.significance) : ""],
        [`${run.name}_direction`, hit ? columnDirection(hit, filters.significance) : ""],
      ];
    })),
  }));
  return <div className="lr-figure">
    <div className="lr-tally">{counts.map((count) => <div key={count.name}>
      <strong>{count.name}</strong>
      <span>{count.opposite} opposite · {count.same} same · {count.total - count.opposite - count.same} no set-level shift</span>
      <small>of {count.total} comparisons, {RULE_LABEL[filters.significance]}</small>
    </div>)}</div>
    <FigureTools view={view} setView={setView} note={`BH and Bonferroni across the engine's ${runs[0].data.multiplicity_family.family_size}-comparison family, before filtering`}
      onCsv={() => csvDownload("comparisons.csv", rows)} />
    {view === "table" ? <div className="lr-table-wrap"><table className="lr-table">
      <thead><tr><th>MoTrPAC comparison</th>{runs.map((run) => <th key={run.name} colSpan={3}>{run.name}</th>)}</tr>
        <tr><th />{runs.map((run) => [<th key={run.name + "m"}>members opposite / measured</th>, <th key={run.name + "t"}>cameraPR t</th>, <th key={run.name + "q"}>{RULE_LABEL[filters.significance].replace(" < 0.05", "")}</th>])}</tr></thead>
      <tbody>{base.map((column) => <tr key={column.column_id}>
        <th scope="row">{column.comparison_label}</th>
        {runs.map((run, index) => {
          const hit = byId[index].get(column.column_id);
          if (!hit) return [<td key={run.name + "a"} colSpan={3}>not measured</td>];
          const label = columnDirection(hit, filters.significance);
          return [
            <td key={run.name + "m"}>{hit.n_opposed}/{hit.n_measured}</td>,
            <td key={run.name + "t"}>{fmt(hit.camera_t)}</td>,
            <td key={run.name + "q"} className={label.startsWith("opposite") ? "lr-pass" : label.startsWith("same") ? "lr-same" : "lr-muted"} title={label}>
              {fmt(columnSignificance(hit, filters.significance))}<small>{label}</small></td>,
          ];
        })}
      </tr>)}</tbody>
    </table></div> : <ColumnDots runs={runs} base={base} byId={byId} rule={filters.significance} />}
  </div>;
}

function ColumnDots({ runs, base, byId, rule }: {
  runs: Array<{ name: string }>; base: ColumnResult[]; byId: Array<Map<string, ColumnResult>>; rule: Significance;
}) {
  const width = 900, top = 56, row = 22, origin = 520, span = 330;
  const height = Math.max(140, top + base.length * row + 20);
  const scale = Math.max(1, ...byId.flatMap((map) => base.map((column) => Math.abs(map.get(column.column_id)?.camera_t ?? 0))));
  const colours = ["#8c2740", "#6d6862"];
  return <div className="lr-chart-scroll"><svg className="lr-plot" width={width} height={height} viewBox={`0 0 ${width} ${height}`} role="img" aria-label="cameraPR t by comparison">
    <rect width={width} height={height} fill="#fbfaf8" />
    <text x="16" y="22" fontSize="13" fontWeight="600" fill="#33312e">cameraPR t: positive = opposite direction to disease</text>
    {runs.map((run, index) => <g key={run.name}><circle cx={16 + index * 200} cy={38} r={5} fill={colours[index]} />
      <text x={26 + index * 200} y={42} fontSize="11" fill="#6d6862">{run.name} (filled: passes {RULE_LABEL[rule]})</text></g>)}
    <line x1={origin} x2={origin} y1={top - 10} y2={height - 10} stroke="#bdb6ae" />
    {base.map((column, i) => {
      const y = top + i * row;
      return <g key={column.column_id}>
        <text x="16" y={y + 4} fontSize="10.5" fill="#33312e">{column.comparison_label}</text>
        {runs.map((run, index) => {
          const hit = byId[index].get(column.column_id);
          if (!hit || hit.camera_t === null) return null;
          const passes = (columnSignificance(hit, rule) ?? 1) < 0.05;
          return <circle key={run.name} cx={origin + (hit.camera_t / scale) * span} cy={y + (index ? 4 : -4)} r={4.5}
            fill={passes ? colours[index] : "#fbfaf8"} stroke={colours[index]} strokeWidth="1.5">
            <title>{`${run.name}: t ${fmt(hit.camera_t)}, ${RULE_LABEL[rule]} value ${fmt(columnSignificance(hit, rule))}`}</title>
          </circle>;
        })}
      </g>;
    })}
  </svg></div>;
}

function useExampleRun(name: string, enabled: boolean): LoadState<AnalysisResponse> {
  const [state, setState] = useState<LoadState<AnalysisResponse>>({ kind: "idle" });
  useEffect(() => {
    if (!enabled || state.kind !== "idle") return;
    setState({ kind: "loading" });
    api.runComparison({ example_name: name, target_species: "human", selected_omics: ["transcriptomics", "proteomics"], fdr_threshold: 0.05 })
      .then((data) => setState({ kind: "ready", data }))
      .catch((error) => setState({ kind: "error", message: errorText(error) }));
  }, [enabled, name, state.kind]);
  return state;
}

function SecondDisease({ filters, active }: { filters: LiveFilters; active: boolean }) {
  const pah = useExampleRun(PAH_EXAMPLE, active);
  const t2d = useExampleRun(T2D_EXAMPLE, active);
  const pending = [pah, t2d].find((state) => state.kind !== "ready");
  return <article className="lr-panel">
    <header className="lr-panel-head">
      <span className="lr-eyebrow">Type 2 diabetes</span>
      <h3>Does a second disease signature behave like the PAH markers in the same MoTrPAC comparisons?</h3>
    </header>
    <div className="lr-contrast-pair">
      <div><span>Rest · disease studies</span><p>PAH: Malenfant 2015, nine proteins lower in PAH muscle. T2D: Mootha 2003, 87 OXPHOS genes lower in impaired-glucose-tolerance or T2D muscle.</p></div>
      <div><span>Exercise · healthy MoTrPAC</span><p>Human muscle after one bout (EE-CON, RE-CON) and rat muscle and heart after 1–8 weeks of training, scored live by the engine.</p></div>
    </div>
    <p className="lr-caveat"><span>Caveat</span>Both signatures are mostly mitochondrial genes lower in disease, so a shared pattern may reflect that gene class rather than either disease. The engine's specificity check (random gene sets of the same class) is in the HTML report.</p>
    {pending?.kind === "error" ? <p className="lr-error">Live engine unavailable: {pending.message}</p>
      : pending ? <p className="lr-small">Running both signatures through the engine…</p>
      : <ColumnTable filters={filters} runs={[
        { name: "PAH", data: (pah as { data: AnalysisResponse }).data },
        { name: "T2D", data: (t2d as { data: AnalysisResponse }).data },
      ]} />}
  </article>;
}

// ---- your data --------------------------------------------------------------------------------

function GeneInput({ filters }: { filters: LiveFilters }) {
  const [text, setText] = useState(GENE_TEMPLATE);
  const [name, setName] = useState("my_signature");
  const [mapping, setMapping] = useState<LoadState<MappingPreview>>({ kind: "idle" });
  const [run, setRun] = useState<LoadState<AnalysisResponse>>({ kind: "idle" });
  const check = () => {
    setMapping({ kind: "loading" });
    api.mappingsPreview({ signature_csv_text: text }).then((data) => setMapping({ kind: "ready", data }))
      .catch((error) => setMapping({ kind: "error", message: errorText(error) }));
  };
  const go = () => {
    setRun({ kind: "loading" });
    api.runComparison({ signature_csv_text: text, signature_name: name, target_species: "human", selected_omics: ["transcriptomics", "proteomics"], fdr_threshold: 0.05 })
      .then((data) => setRun({ kind: "ready", data }))
      .catch((error) => setRun({ kind: "error", message: errorText(error) }));
  };
  const unmapped = mapping.kind === "ready" ? mapping.data.rows.filter((row) => row.status !== "mapped") : [];
  return <article className="lr-panel">
    <header className="lr-panel-head"><span className="lr-eyebrow">Gene or protein set</span>
      <h3>Score your own disease signature against every MoTrPAC muscle comparison</h3></header>
    <p className="lr-small">One row per gene: <code>gene_symbol,direction</code> with +1 higher or −1 lower in disease. UniProt, Ensembl and rat symbols are also accepted. The engine covers human muscle (endurance, resistance × 3 times) and rat gastrocnemius, vastus lateralis and heart (1–8 weeks, by sex).</p>
    <div className="lr-input-row">
      <label>Name<input value={name} onChange={(event) => setName(event.target.value)} /></label>
      <label className="lr-file">Upload CSV<input type="file" accept=".csv,.txt" onChange={(event) => readFile(event.target.files?.[0], setText)} /></label>
    </div>
    <textarea className="lr-textarea" value={text} onChange={(event) => { setText(event.target.value); setMapping({ kind: "idle" }); }} rows={7} spellCheck={false} />
    <div className="lr-actions">
      <button type="button" className="lr-outline-button" onClick={check} disabled={mapping.kind === "loading"}>Check IDs</button>
      <button type="button" className="lr-primary-button" onClick={go} disabled={run.kind === "loading"}>{run.kind === "loading" ? "Running…" : "Run comparison"}</button>
    </div>
    {mapping.kind === "error" && <p className="lr-error">{mapping.message}</p>}
    {mapping.kind === "ready" && <p className="lr-small">
      {mapping.data.counts.n_mapped ?? 0} mapped · {mapping.data.counts.n_ambiguous ?? 0} ambiguous · {mapping.data.counts.n_unmapped ?? 0} not found
      {unmapped.length > 0 && <> · <strong>Not mapped:</strong> {unmapped.map((row) => row.input_id).join(", ")}</>}
    </p>}
    {run.kind === "error" && <p className="lr-error">{run.message}</p>}
    {run.kind === "ready" && <>
      <p className="lr-small">{run.data.n_counted_genes} of {run.data.n_input_rows} rows counted · run <code>{run.data.run_id}</code> ·{" "}
        <a href={api.reportUrl(run.data.run_id)} target="_blank" rel="noreferrer">HTML report</a> · <a href={api.exportUrl(run.data.run_id)}>evidence bundle</a></p>
      <ColumnTable filters={filters} runs={[{ name: run.data.signature_name, data: run.data }]} />
    </>}
  </article>;
}

function MetaboliteInput() {
  const [text, setText] = useState(METABOLITE_TEMPLATE);
  const [state, setState] = useState<LoadState<MetabolitePreview>>({ kind: "idle" });
  const go = () => {
    setState({ kind: "loading" });
    api.metabolitePreview({ csv_text: text }).then((data) => setState({ kind: "ready", data }))
      .catch((error) => setState({ kind: "error", message: errorText(error) }));
  };
  return <article className="lr-panel">
    <header className="lr-panel-head"><span className="lr-eyebrow">Metabolite list</span>
      <h3>Check whether reported disease metabolites are stable in healthy blood</h3></header>
    <p className="lr-small">Columns as in <code>MoTrPAC Hackathon/Metabolomics/README.md</code>. Names match MoTrPAC by exact RefMet name. A context label needs BH q &lt; 0.05 and |log2 FC| ≥ 0.5.</p>
    <div className="lr-input-row">
      <label className="lr-file">Upload CSV<input type="file" accept=".csv,.txt" onChange={(event) => readFile(event.target.files?.[0], setText)} /></label>
      <button type="button" className="lr-text-button" onClick={() => download("metabolite_list_template.csv", new Blob([METABOLITE_TEMPLATE], { type: "text/csv" }))}>Download template</button>
    </div>
    <textarea className="lr-textarea" value={text} onChange={(event) => setText(event.target.value)} rows={5} spellCheck={false} />
    <div className="lr-actions"><button type="button" className="lr-primary-button" onClick={go} disabled={state.kind === "loading"}>Match to MoTrPAC</button></div>
    {state.kind === "error" && <p className="lr-error">{state.message}</p>}
    {state.kind === "ready" && <>
      <p className="lr-small">{state.data.n_blood_matched} of {state.data.n_input_rows} matched in MoTrPAC blood
        {state.data.unmapped_ids.length > 0 && <> · <strong>Not matched:</strong> {state.data.unmapped_ids.join(", ")}</>}</p>
      <div className="lr-table-wrap"><table className="lr-table">
        <thead><tr><th>Metabolite</th><th>Disease direction</th><th>Context</th><th>Blood / muscle match</th><th>MoTrPAC rows</th></tr></thead>
        <tbody>{state.data.rows.map((row) => <tr key={row.input_row}>
          <th scope="row">{row.name}<small>{row.list_id}</small></th><td>{row.pah_direction}</td><td>{row.context_label}</td>
          <td>{row.blood_matched ? "yes" : "no"} / {row.muscle_matched ? "yes" : "no"}</td><td>{row.observations.length}</td>
        </tr>)}</tbody></table></div>
      <div className="lr-actions"><button type="button" className="lr-outline-button" onClick={() => csvDownload("metabolite_context.csv",
        state.data.rows.flatMap((row) => row.observations.map((obs) => ({ name: row.name, context: row.context_label, ...obs }))))}>CSV of all MoTrPAC rows</button></div>
      <p className="lr-small">{state.data.note}</p>
    </>}
  </article>;
}

// ---- catalog ------------------------------------------------------------------------------------

const CATALOG_PANELS = ["human_muscle_protein", "human_muscle_rna", "blood_go", "blood_biocarta", "rat_muscle_protein", "metabolomics_context"];

function Catalog({ story, filters }: { story: StoryResponse; filters: LiveFilters }) {
  const [limit, setLimit] = useState(60);
  const records = useMemo(() => filterRecords(story.panels
    .filter((panel) => CATALOG_PANELS.includes(panel.id))
    .flatMap((panel) => panel.records)
    .filter((record) => record.molecule && !record.record_type.includes("summary")), filters), [story, filters]);
  const counts = Object.values(AGREEMENT).map((label) => [label, records.filter((record) => agreementClass(record, filters.significance) === label).length] as const);
  return <article className="lr-panel">
    <header className="lr-panel-head"><span className="lr-eyebrow">Evidence catalog</span>
      <h3>Every molecule, layer, tissue and time in one table</h3></header>
    <p className="lr-small">Direction of the disease change versus the healthy MoTrPAC change, under {RULE_LABEL[filters.significance]}. Rows are cross-study comparisons, not treatment tests. This table is the input for a later discordance classifier.</p>
    <div className="lr-tally">{counts.map(([label, count]) => <div key={label}><strong>{count}</strong><span>{label}</span></div>)}</div>
    <div className="lr-figure-tools"><span>{records.length} rows</span>
      <button type="button" className="lr-outline-button" onClick={() => csvDownload("evidence_catalog.csv", records.map((record) => recordRow(record, filters.significance)))}>CSV</button></div>
    <div className="lr-table-wrap"><table className="lr-table">
      <thead><tr><th>Molecule</th><th>Context</th><th>Disease</th><th>MoTrPAC</th><th>Estimate</th><th>{RULE_LABEL[filters.significance].replace(" < 0.05", "")}</th><th>Class</th></tr></thead>
      <tbody>{records.slice(0, limit).map((record) => {
        const value = significanceValue(record, filters.significance);
        return <tr key={record.id}>
          <th scope="row">{record.molecule}<small>{record.panel_id.replace(/_/g, " ")}</small></th>
          <td>{record.species} · {record.tissue} · {record.layer}<small>{optionLabel(record.contrast ?? "")} · {optionLabel(record.timepoint ?? "")}{record.sex && !["all", "both"].includes(record.sex) ? " · " + record.sex : ""}</small></td>
          <td>{record.disease_direction ?? "—"}</td><td>{record.motrpac_measured === false ? "not measured" : record.motrpac_direction ?? "—"}</td>
          <td>{fmt(record.value)}</td>
          <td title={value === null ? unavailableRuleReason(record, filters.significance) : undefined}>{value === null ? "—" : fmt(value)}</td>
          <td>{agreementClass(record, filters.significance)}</td>
        </tr>;
      })}</tbody>
    </table></div>
    {records.length > limit && <div className="lr-actions"><button type="button" className="lr-outline-button" onClick={() => setLimit(limit + 200)}>Show 200 more</button></div>}
  </article>;
}

// ---- section ------------------------------------------------------------------------------------

export default function LiveResultsStory() {
  const [story, setStory] = useState<LoadState<StoryResponse>>({ kind: "loading" });
  const [filters, setFilters] = useState<LiveFilters>(readFilters);
  const [step, setStep] = useState<StepId>(readStep);

  const load = () => {
    setStory({ kind: "loading" });
    api.story().then((data) => setStory({ kind: "ready", data }))
      .catch((error) => setStory({ kind: "error", message: errorText(error) }));
  };
  useEffect(load, []);
  useEffect(() => writeUrl(filters, step), [filters, step]);

  const contexts = useMemo<FilterContext[]>(() => story.kind === "ready"
    ? story.data.panels.flatMap((panel) => panel.records.map((record) => ({
      species: record.species, tissue: record.tissue, layer: record.layer, exercise: record.exercise,
      contrast: record.contrast, timepoint: record.timepoint, sex: record.sex, motrpac_measured: record.motrpac_measured,
    }))) : [], [story]);

  const showPanelData = (context: FilterContext) => setFilters({ ...DEFAULT_FILTERS, significance: filters.significance,
    species: context.species ?? ALL });

  return <section className="view lr-view">
    <div className="lr-intro">
      <span className="lr-intro-number">05</span>
      <div>
        <span className="lr-eyebrow">Live results</span>
        <h2>PAH markers in healthy exercise physiology</h2>
        <p>Published PAH results, then the same molecules in MoTrPAC: protein, RNA, rat training, metabolites, and a second disease. Every number comes from a committed output or a live engine run.</p>
      </div>
    </div>
    <p className="lr-disclaimer"><span>Separate cohorts</span>MoTrPAC people and rats are healthy. A same or opposite direction is a question for a patient study, not a treatment effect.</p>

    {story.kind === "loading" && <p className="lr-small">Loading results…</p>}
    {story.kind === "error" && <div className="lr-empty"><strong>The results API is not reachable.</strong>
      <p>{story.message}</p><p>Start it as described in <code>RUN_LOCAL.md</code>, then retry.</p>
      <button type="button" className="lr-outline-button" onClick={load}>Retry</button></div>}

    {story.kind === "ready" && <>
      <nav className="lr-steps" aria-label="Story steps">
        {STEPS.map((item, index) => <button type="button" key={item.id} aria-current={step === item.id ? "step" : undefined}
          className={"lr-step" + (step === item.id ? " lr-step--active" : "")} onClick={() => setStep(item.id)}>
          <span>{index < 5 ? index + 1 : item.id === "data" ? "+" : "≡"}</span>
          <strong>{item.title}</strong><small>{item.subtitle}</small>
        </button>)}
      </nav>
      <FilterBar contexts={contexts} filters={filters} onChange={setFilters} onReset={() => setFilters({ ...DEFAULT_FILTERS, significance: filters.significance })} />

      <div className="lr-step-body">
        {step === "2" && <LayerDiscordance story={story.data} filters={filters} />}
        {STEPS.find((item) => item.id === step)!.panels.map((id) => {
          const panel = story.data.panels.find((item) => item.id === id);
          return panel ? <PanelCard key={id} panel={panel} filters={filters} onUseContext={showPanelData} /> : null;
        })}
        {step === "5" && <SecondDisease filters={filters} active={step === "5"} />}
        {step === "data" && <><GeneInput filters={filters} /><MetaboliteInput /></>}
        {step === "catalog" && <><LayerDiscordance story={story.data} filters={filters} /><Catalog story={story.data} filters={filters} /></>}
      </div>

      <nav className="lr-pager" aria-label="Previous and next step">
        {STEPS.findIndex((item) => item.id === step) > 0 && <button type="button" className="lr-outline-button"
          onClick={() => setStep(STEPS[STEPS.findIndex((item) => item.id === step) - 1].id)}>← Previous</button>}
        <span />
        {STEPS.findIndex((item) => item.id === step) < STEPS.length - 1 && <button type="button" className="lr-primary-button"
          onClick={() => { setStep(STEPS[STEPS.findIndex((item) => item.id === step) + 1].id); window.scrollTo({ top: 0, behavior: "smooth" }); }}>
          Next: {STEPS[STEPS.findIndex((item) => item.id === step) + 1].title} →</button>}
      </nav>
      <p className="lr-footnote">MoTrPAC human data: MotrpacHumanPreSuspensionAnalysis {story.data.provenance?.human_package_version} ({story.data.provenance?.human_collection}). Sources are listed under each figure.</p>
    </>}
  </section>;
}
