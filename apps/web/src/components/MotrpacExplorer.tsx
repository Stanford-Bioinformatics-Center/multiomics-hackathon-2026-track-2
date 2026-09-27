import { useEffect, useMemo, useRef, useState, type ReactNode, type RefObject } from "react";
import {
  api, ApiError, type ExplorerColumn, type ExplorerExample, type ExplorerKind, type ExplorerLayer,
  type ExplorerResponse, type ExplorerValue,
} from "../api/client";
import { Bars, COLOURS, DotPlot, HeatGrid, Scatter, SmallMultiples, type Cell, type Series } from "./liveCharts";
import { Pick, RULE_LABEL, RuleSwitch, csvDownload, fmt, svgPng, type Significance } from "./ui";
import "./LiveDashboard.css";

type Load<T> = { kind: "idle" } | { kind: "loading" } | { kind: "error"; message: string } | { kind: "ready"; data: T };
type Filters = { species: string; tissue: string; contrast: string; sex: string; time: string; layer: string };

const KINDS: Array<{ id: ExplorerKind; label: string; hint: string; placeholder: string }> = [
  { id: "genes", label: "Genes & proteins", hint: "One per line, or columns gene_symbol / uniprot / ensembl. Optional direction (+1, −1) or score.",
    placeholder: "gene_symbol,direction\nNDUFA9,-1\nUQCRC1,-1\nPPARGC1A,1" },
  { id: "metabolites", label: "Metabolites", hint: "RefMet names (or hmdb / kegg ids). Optional direction or log2fc.",
    placeholder: "refmet_name,log2fc\nSuccinic acid,0.6\nLactic acid,0.4\nOleic acid,-1.1" },
  { id: "pathways", label: "Pathways", hint: "Pathway names, optionally with a source column (BioCarta, KEGG, Reactome, GO…).",
    placeholder: "pathway,source\nTCR PATHWAY,BioCarta\nOXIDATIVE PHOSPHORYLATION,GOBP" },
];
const CATEGORY: Record<string, string> = {
  "EE-CON": "Endurance vs control", "RE-CON": "Resistance vs control", "EE-EE": "Endurance, before → after",
  "RE-RE": "Resistance, before → after", "CON-CON": "Control, no exercise", "TRAIN-SED": "Trained vs sedentary",
  "VS-CON": "Exercise vs control (both)",
};
const CATEGORY_ORDER = ["EE-CON", "RE-CON", "EE-EE", "RE-RE", "CON-CON", "TRAIN-SED"];
const TISSUE_ORDER = ["VL", "muscle", "BLOOD", "blood", "ADIPOSE", "adipose", "SKM-GN", "SKM-VL", "HEART", "WAT-SC", "BAT", "PLASMA", "LIVER", "LUNG", "KIDNEY"];
const ALL = "all";
const SERIES_COLOURS = [COLOURS.warm, COLOURS.cool, "#1f8a6b", "#b07a16", "#6d6862"];

const rank = (list: string[], v: string) => { const i = list.indexOf(v); return i === -1 ? 999 : i; };
const uniq = <T,>(xs: T[]) => Array.from(new Set(xs));
const shortTime = (t: string) => t.replace(" after", "").replace(" (during)", " during");
function errorText(error: unknown): string { return error instanceof ApiError ? `${error.status}: ${error.detail}` : String(error); }
function setP(c: ExplorerColumn, rule: Significance, ranked: boolean): number | null | undefined {
  if (ranked) return rule === "raw" ? c.rho_p : rule === "bh" ? c.rho_bh : c.rho_bonferroni;
  return rule === "raw" ? c.set_p : rule === "bh" ? c.set_bh : c.set_bonferroni;
}
const valueP = (v: ExplorerValue, rule: Significance) => rule === "raw" ? v[3] : rule === "bh" ? v[4] : v[5];
const passes = (p: number | null | undefined) => p === null || p === undefined ? null : p < 0.05;
const inContrast = (c: ExplorerColumn, contrast: string) => contrast === ALL || c.category === contrast || (contrast === "VS-CON" && (c.category === "EE-CON" || c.category === "RE-CON"));

async function toBase64(file: File): Promise<string> {
  const bytes = new Uint8Array(await file.arrayBuffer());
  let binary = "";
  for (let i = 0; i < bytes.length; i += 0x8000) binary += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
  return btoa(binary);
}

// ---- input panel -------------------------------------------------------------------------------------

type Draft = { text: string; example: ExplorerExample | null; name: string | null; file: string | null };
const EMPTY: Draft = { text: "", example: null, name: null, file: null };

function InputPanel({ examples, drafts, setDrafts, onRun, running }: {
  examples: ExplorerExample[]; drafts: Record<ExplorerKind, Draft>; setDrafts: (d: Record<ExplorerKind, Draft>) => void;
  onRun: () => void; running: boolean;
}) {
  const [kind, setKind] = useState<ExplorerKind>("genes");
  const [uploadError, setUploadError] = useState<string | null>(null);
  const meta = KINDS.find((k) => k.id === kind)!;
  const draft = drafts[kind];
  const update = (next: Draft) => setDrafts({ ...drafts, [kind]: next });
  const upload = async (file: File | undefined) => {
    if (!file) return;
    setUploadError(null);
    try {
      const read = await api.explorerReadFile({ filename: file.name, content_base64: await toBase64(file) });
      const target = read.kind ?? kind;
      setKind(target);
      // A new upload replaces any previous list, so each analysis shows one input.
      setDrafts({ genes: EMPTY, metabolites: EMPTY, pathways: EMPTY, [target]: { text: read.text, example: null, name: file.name.replace(/\.[^.]+$/, ""), file: `${file.name} · ${read.n_rows} rows` } });
    } catch (e) { setUploadError(errorText(e)); }
  };
  const ready = KINDS.filter((k) => drafts[k.id].example || drafts[k.id].text.trim());
  return <article className="lr-card xp-input">
    <div className="xp-input-top">
      <div className="lr-switch" role="tablist" aria-label="List type">
        {KINDS.map((k) => <button type="button" key={k.id} aria-pressed={kind === k.id} onClick={() => setKind(k.id)}>
          {k.label}{(drafts[k.id].example || drafts[k.id].text.trim()) && <span className="xp-dot" aria-label="has a list" />}</button>)}
      </div>
      <label className="lr-button xp-upload">Upload spreadsheet<input type="file" accept=".xlsx,.xls,.csv,.tsv,.txt"
        onChange={(e) => { upload(e.target.files?.[0]); e.currentTarget.value = ""; }} /></label>
    </div>
    <div className="xp-input-grid">
      <div>
        {draft.example
          ? <div className="xp-example-chosen"><strong>{draft.example.label}</strong><span>{draft.example.source}</span>
              <button type="button" className="lr-link" onClick={() => update(EMPTY)}>Use my own list instead</button></div>
          : <>
              {draft.file && <div className="xp-file"><strong>{draft.file}</strong><button type="button" className="lr-link" onClick={() => update(EMPTY)}>Clear</button></div>}
              <textarea className="lr-textarea" value={draft.text} placeholder={meta.placeholder} spellCheck={false}
                aria-label={meta.label} onChange={(e) => setDrafts({ genes: EMPTY, metabolites: EMPTY, pathways: EMPTY, [kind]: { ...draft, text: e.target.value, example: null } })} rows={8} />
            </>}
        <p className="lr-hint">{meta.hint}</p>
        {uploadError && <p className="lr-error">{uploadError}</p>}
      </div>
      <div className="xp-examples">
        <span className="xp-label">Examples</span>
        {examples.filter((e) => e.kind === kind).map((e) => <button type="button" key={e.id}
          className={"xp-example" + (draft.example?.id === e.id ? " is-on" : "")} onClick={() => setDrafts({ genes: EMPTY, metabolites: EMPTY, pathways: EMPTY, [kind]: { ...EMPTY, example: e } })}>
          <strong>{e.label}</strong><small>{e.source}</small></button>)}
      </div>
    </div>
    <footer className="xp-input-foot">
      <span className="lr-hint">{ready.length ? `Ready: ${ready.map((k) => k.label.toLowerCase()).join(", ")}` : "Upload a spreadsheet, paste a list, or pick an example."}</span>
      <button type="button" className="lr-button" disabled={!ready.length || running} onClick={onRun}>{running ? "Matching to MoTrPAC…" : "Find in MoTrPAC"}</button>
    </footer>
  </article>;
}

// ---- shared card shell ---------------------------------------------------------------------------------

function Card({ title, chips, stat, caption, children, foot }: {
  title: string; chips: string[]; stat?: string; caption?: string; children: ReactNode; foot?: ReactNode;
}) {
  return <article className="lr-card">
    <header className="lr-card-head">
      <div className="lr-card-text"><h3>{title}</h3><div className="lr-chips">{chips.map((c) => <span key={c}>{c}</span>)}</div></div>
      {stat && <div className="lr-stat"><strong>{stat}</strong><span>{caption}</span></div>}
    </header>
    <div className="lr-visual">{children}</div>
    {foot && <footer className="lr-card-foot">{foot}</footer>}
  </article>;
}
const PngButton = ({ svgRef, name }: { svgRef: RefObject<SVGSVGElement | null>; name: string }) =>
  <button type="button" className="lr-ghost" onClick={() => svgPng(svgRef.current, name).catch(() => undefined)}>PNG</button>;

// ---- gene / protein / metabolite layers ---------------------------------------------------------------

const TISSUE_HINTS: Array<[RegExp, string[]]> = [
  [/plasma|serum|blood|pbmc/i, ["BLOOD", "PLASMA", "blood"]],
  [/muscle|vastus|gastroc/i, ["VL", "SKM-GN", "SKM-VL", "muscle"]],
  [/adipose|fat/i, ["ADIPOSE", "WAT-SC", "adipose"]],
];

/** Default view: human if measured; the tissue named in the list's name, else the one measuring most of its molecules. */
function defaults(cols: ExplorerColumn[], values: ExplorerValue[], allCols: ExplorerColumn[], species?: string, hint?: string): Filters {
  const sp = species ?? (cols.some((c) => c.species === "human") ? "human" : "rat");
  const seen = new Map<string, Set<number>>();
  for (const v of values) {
    const c = allCols[v[1]];
    if (c && c.species === sp) { if (!seen.has(c.tissue)) seen.set(c.tissue, new Set()); seen.get(c.tissue)!.add(v[0]); }
  }
  const coverage = (t: string) => seen.get(t)?.size ?? 0;
  const tissues = uniq(cols.filter((c) => c.species === sp).map((c) => c.tissue)).sort((a, b) => coverage(b) - coverage(a) || rank(TISSUE_ORDER, a) - rank(TISSUE_ORDER, b));
  const cats = uniq(cols.filter((c) => c.species === sp).map((c) => c.category));
  const hinted = hint ? TISSUE_HINTS.find(([re]) => re.test(hint))?.[1].find((t) => tissues.includes(t)) : undefined;
  return { species: sp, tissue: hinted ?? tissues[0] ?? ALL, contrast: cats.includes("EE-CON") ? "EE-CON" : cats[0] ?? ALL, sex: ALL, time: ALL, layer: ALL };
}

function LayerView({ layer, directed, ranked, rule, hint }: { layer: ExplorerLayer; directed: boolean; ranked: boolean; rule: Significance; hint?: string }) {
  const cols = layer.columns ?? [];
  const values = layer.values ?? [];
  const molecules = layer.molecules ?? [];
  const [f, setF] = useState<Filters>(() => defaults(cols, values, cols, undefined, hint));
  useEffect(() => setF(defaults(cols, values, cols, undefined, hint)), [layer]); // eslint-disable-line react-hooks/exhaustive-deps
  const refOverview = useRef<SVGSVGElement>(null), refMolecules = useRef<SVGSVGElement>(null), refTime = useRef<SVGSVGElement>(null);
  const refContext = useRef<SVGSVGElement>(null), refBars = useRef<SVGSVGElement>(null);

  const colIndex = useMemo(() => new Map(cols.map((c, j) => [c, j])), [cols]);
  const valueAt = useMemo(() => { const m = new Map<string, ExplorerValue>(); for (const v of values) m.set(`${v[0]}|${v[1]}`, v); return m; }, [values]);
  const human = f.species === "human";
  const isMetab = layer.id === "metabolomics";
  const bySpecies = cols.filter((c) => c.species === f.species);
  const tissues = uniq(bySpecies.map((c) => c.tissue)).sort((a, b) => rank(TISSUE_ORDER, a) - rank(TISSUE_ORDER, b));
  const tissueLabel = (t: string) => cols.find((c) => c.tissue === t)?.tissue_label ?? t;
  const inTissue = bySpecies.filter((c) => c.tissue === f.tissue);
  const categories = uniq(inTissue.map((c) => c.category));
  const selected = inTissue.filter((c) => inContrast(c, f.contrast) && (f.sex === ALL || c.sex === f.sex) && (f.time === ALL || c.time === f.time));
  const timesOf = (list: ExplorerColumn[]) => uniq(list.map((c) => c.time)).sort((a, b) => cols.find((c) => c.time === a)!.time_rank - cols.find((c) => c.time === b)!.time_rank);
  const timeLabel = (t: string) => cols.find((c) => c.time === t)?.time_label ?? t;
  const rowName = (i: number) => `${molecules[i].name}${molecules[i].dir > 0 ? "  ↑" : molecules[i].dir < 0 ? "  ↓" : ""}`;
  const setValue = (c: ExplorerColumn) => ranked ? c.rho ?? null : c.set_t ?? null;
  const setLabels: [string, string] = ranked ? ["opposite ranking", "same ranking"] : directed ? ["same direction as list", "opposite"] : ["down", "up"];
  const describe = (c: ExplorerColumn) => `${c.tissue_label} · ${CATEGORY[c.category]} · ${c.time_label}${c.sex !== "all" ? " · " + c.sex : ""}`;
  const rows = molecules.map((_, i) => i);

  // overview: every comparison for this species
  const multiLayer = uniq(bySpecies.map((c) => c.layer)).length > 1;
  const overviewRows = uniq(bySpecies.map((c) => human ? `${c.tissue}|${c.category}|${c.layer}` : `${c.tissue}|${c.layer}`)).sort((a, b) => {
    const [ta, ca] = a.split("|"), [tb, cb] = b.split("|");
    return rank(TISSUE_ORDER, ta) - rank(TISSUE_ORDER, tb) || rank(CATEGORY_ORDER, ca) - rank(CATEGORY_ORDER, cb) || a.localeCompare(b);
  });
  const overviewLabel = (key: string) => { const [t, c, l] = key.split("|"); return human ? `${tissueLabel(t)} · ${CATEGORY[c]}${multiLayer ? ` (${l})` : ""}` : `${tissueLabel(t)}${multiLayer ? ` (${c})` : ""}`; };
  const allTimes = timesOf(bySpecies);
  const overviewGroups = human ? [{ cols: allTimes.map((t) => ({ key: t, label: shortTime(timeLabel(t)) })) }]
    : ["female", "male"].map((s) => ({ title: s === "female" ? "Female" : "Male", cols: allTimes.map((t) => ({ key: `${s}|${t}`, label: t })) }));
  const overviewFind = (rowKey: string, colKey: string) => {
    const [t, c, l] = rowKey.split("|");
    if (human) return bySpecies.find((x) => x.tissue === t && x.category === c && x.layer === l && x.time === colKey);
    const [s, tm] = colKey.split("|");
    return bySpecies.find((x) => x.tissue === t && x.layer === c && x.sex === s && x.time === tm);
  };
  const tested = bySpecies.filter((c) => setP(c, rule, ranked) !== null && setP(c, rule, ranked) !== undefined);
  const overviewSig = tested.filter((c) => passes(setP(c, rule, ranked))).length;

  // molecules x time for the selection
  const groupKey = (c: ExplorerColumn) => !human && f.sex === ALL ? (c.sex === "female" ? "Female" : "Male") : f.contrast === ALL || f.contrast === "VS-CON" ? CATEGORY[c.category] : "";
  const moleculeGroups = uniq(selected.map(groupKey)).map((title) => ({ title: title || undefined,
    cols: selected.filter((c) => groupKey(c) === title).sort((a, b) => a.time_rank - b.time_rank).map((c) => ({ key: String(colIndex.get(c)), label: shortTime(c.time_label) })) }));
  const moleculeCell = (row: string, key: string): Cell => {
    const i = rows.find((k) => rowName(k) === row)!;
    const v = valueAt.get(`${i}|${key}`);
    if (!v) return null;
    return { value: v[2], pass: passes(valueP(v, rule)), title: `${molecules[i].name} · ${describe(cols[Number(key)])}\nlog2 FC ${fmt(v[2])} · raw P ${fmt(v[3])} · BH q ${fmt(v[4])} · Bonferroni ${fmt(v[5])}` };
  };
  const changed = rows.filter((i) => selected.some((c) => passes(valueP(valueAt.get(`${i}|${colIndex.get(c)}`) ?? [0, 0, null, null, null, null], rule)))).length;

  // over time: one panel per molecule, one line per sex (rat) or contrast (human)
  const lineKey = (c: ExplorerColumn) => human ? CATEGORY[c.category] : c.sex === "female" ? "Female" : "Male";
  const lineNames = uniq(selected.map(lineKey));
  const lineTimes = timesOf(selected);
  const timePanels = rows.slice(0, 12).map((i) => ({
    title: molecules[i].name,
    series: lineNames.map((name, k): Series => ({ name, colour: SERIES_COLOURS[k % SERIES_COLOURS.length], points: lineTimes.map((t, x) => {
      const c = selected.find((col) => lineKey(col) === name && col.time === t);
      const v = c ? valueAt.get(`${i}|${colIndex.get(c)}`) : undefined;
      return { x, y: v ? v[2] : null, pass: v ? passes(valueP(v, rule)) : null,
        title: c && v ? `${molecules[i].name} · ${describe(c)}\nlog2 FC ${fmt(v[2])} · ${RULE_LABEL[rule].split(" <")[0]} ${fmt(valueP(v, rule))}` : "" };
    }) })),
  }));

  // metabolites in healthy people: exercise effect vs drift in resting controls
  const context = useMemo(() => {
    if (!isMetab || !human) return null;
    const tissueCols = bySpecies.filter((c) => c.tissue === f.tissue);
    const ex = tissueCols.filter((c) => c.category === "EE-CON"), drift = tissueCols.filter((c) => c.category === "CON-CON");
    const big = (v: ExplorerValue) => passes(valueP(v, rule)) === true && Math.abs(v[2] ?? 0) >= 0.5;
    const out = rows.map((i) => {
      const exV = ex.map((c) => valueAt.get(`${i}|${colIndex.get(c)}`)).filter(Boolean) as ExplorerValue[];
      const drV = drift.map((c) => valueAt.get(`${i}|${colIndex.get(c)}`)).filter(Boolean) as ExplorerValue[];
      if (!exV.length && !drV.length) return { i, label: "not measured here", drift: null as number | null };
      const label = exV.some(big) ? "exercise-sensitive" : drV.some(big) ? "drifts without exercise" : "stable";
      return { i, label, drift: drV.length ? Math.max(...drV.map((v) => Math.abs(v[2] ?? 0))) : null };
    });
    const withScore = out.filter((o) => o.drift !== null && o.drift > 0 && molecules[o.i].score !== null);
    const ratios = withScore.map((o) => Math.abs(molecules[o.i].score!) / o.drift!).sort((a, b) => a - b);
    return { out, withScore, median: ratios.length ? ratios[Math.floor(ratios.length / 2)] : null };
  }, [isMetab, human, f.tissue, rule, layer]); // eslint-disable-line react-hooks/exhaustive-deps

  const barItems = selected.filter((c) => setValue(c) !== null).sort((a, b) => a.time_rank - b.time_rank || a.sex.localeCompare(b.sex)).slice(0, 16).map((c) => ({
    label: `${shortTime(c.time_label)}${!human && f.sex === ALL ? (c.sex === "female" ? " F" : " M") : ""}${(f.contrast === ALL || f.contrast === "VS-CON") && human ? " " + c.category : ""}`,
    value: setValue(c), pass: passes(setP(c, rule, ranked)),
    title: `${describe(c)}\n${ranked ? "Spearman ρ" : "cameraPR t"} ${fmt(setValue(c))} · ${RULE_LABEL[rule].split(" <")[0]} ${fmt(setP(c, rule, ranked))}`,
  }));
  const pick = (key: keyof Filters, value: string) => setF(key === "species" ? defaults(cols.filter((c) => c.species === value), values, cols, value, hint)
    : { ...f, [key]: value, ...(key === "tissue" ? { time: ALL } : {}) });
  const large = layer.n_matched > (layer.n_shown ?? 0);
  const listSize = `${layer.n_shown}${large ? ` of ${layer.n_matched.toLocaleString()} (${ranked ? "largest scores" : "first in your list"})` : ""} molecules`;
  const cameraNote = ranked ? "Spearman ρ across shared genes" : directed ? "cameraPR t · positive = opposite to your list" : "cameraPR t · positive = up";

  const cards: Record<string, ReactNode> = {
    overview: (overviewRows.length > 0 && <Card title="Every MoTrPAC comparison"
      chips={[human ? "Human · one exercise bout" : "Rat · endurance training", ranked ? "Rank correlation with your scores" : "Set-level test of your list (cameraPR)", "Click a cell to open it"]}
      stat={`${overviewSig}/${tested.length}`} caption={`comparisons pass ${RULE_LABEL[rule]}`}
      foot={<><span className="lr-hint">One multiple-testing family: {layer.family?.scope}{layer.family?.n_tests ? ` · ${layer.family.n_tests} tests` : ""}.{ranked ? " Rank-correlation P values treat genes as independent, so read ρ as descriptive." : ""}</span><PngButton svgRef={refOverview} name={`${layer.id}_overview.png`} /></>}>
      <HeatGrid svgRef={refOverview} rows={overviewRows.map(overviewLabel)} rowLabelWidth={300} groups={overviewGroups} labels={setLabels}
        cell={(label, key) => { const c = overviewFind(overviewRows.find((r) => overviewLabel(r) === label)!, key);
          return c ? { value: setValue(c), pass: passes(setP(c, rule, ranked)), title: `${describe(c)}\n${ranked ? "ρ" : "t"} ${fmt(setValue(c))} · ${RULE_LABEL[rule].split(" <")[0]} ${fmt(setP(c, rule, ranked))}` } : null; }}
        onCell={(label, key) => { const c = overviewFind(overviewRows.find((r) => overviewLabel(r) === label)!, key);
          if (c) setF({ ...f, species: c.species, tissue: c.tissue, contrast: c.category, sex: human ? ALL : c.sex, time: ALL }); }} />
    </Card>),
    molecules: (<Card title={`Your ${isMetab ? "metabolites" : "molecules"} in ${tissueLabel(f.tissue).toLowerCase()}`}
      chips={[f.contrast === ALL ? "All contrasts" : CATEGORY[f.contrast], listSize, ...(directed ? ["↑ ↓ = direction in your list"] : [])]}
      stat={selected.length ? `${changed}/${rows.length}` : undefined} caption={`change significantly at ≥ 1 time (${RULE_LABEL[rule]})`}
      foot={<><span className="lr-hint">{layer.molecule_rule}</span>
        <button type="button" className="lr-ghost" onClick={() => csvDownload(`${layer.id}_molecules.csv`, values.map((v) => {
          const c = cols[v[1]];
          return { molecule: molecules[v[0]].name, list_direction: molecules[v[0]].dir, species: c.species, tissue: c.tissue_label, layer: c.layer, contrast: c.category, time: c.time_label, sex: c.sex, log2fc: v[2], raw_p: v[3], bh_q: v[4], bonferroni: v[5] };
        }))}>CSV</button><PngButton svgRef={refMolecules} name={`${layer.id}_molecules.png`} /></>}>
      {selected.length && moleculeGroups.length ? <HeatGrid svgRef={refMolecules} rows={rows.map(rowName)} groups={moleculeGroups} rowLabelWidth={isMetab ? 220 : 170}
        cell={moleculeCell} minScale={0.5} labels={["lower after exercise", "higher"]} /> : <div className="lr-empty">No MoTrPAC comparison matches these filters.</div>}
    </Card>),
    time: (!ranked && timePanels.length > 0 && lineTimes.length > 1 && <Card title="Over time"
      chips={[tissueLabel(f.tissue), human ? (f.contrast === ALL ? "All contrasts" : CATEGORY[f.contrast]) : "Trained vs sedentary, by sex", rows.length > 12 ? "First 12 molecules" : `${rows.length} molecules`]}
      foot={<PngButton svgRef={refTime} name={`${layer.id}_over_time.png`} />}>
      <SmallMultiples svgRef={refTime} panels={timePanels} xLabels={lineTimes.map((t) => shortTime(timeLabel(t)))} yLabel={human ? "log2 FC vs control" : "log2 FC trained vs sedentary"} />
    </Card>),
    context: (context && <Card title={`How these metabolites behave in healthy ${tissueLabel(f.tissue).toLowerCase()}`}
      chips={["Exercise effect: endurance vs control", "Drift: controls, no exercise", "A change needs significance and |log2 FC| ≥ 0.5"]}
      stat={context.median !== null ? `${context.median.toFixed(1)}×` : undefined} caption="median: your difference vs largest healthy drift"
      foot={<PngButton svgRef={refContext} name={`${layer.id}_context.png`} />}>
      <div className="xp-context">
        <div className="lr-tiles">{["stable", "exercise-sensitive", "drifts without exercise", "not measured here"].map((label) => {
          const n = context.out.filter((o) => o.label === label).length;
          return n ? <div className="lr-tile" key={label}><strong>{n}</strong><span>{label}</span></div> : null;
        })}</div>
        {context.withScore.length > 0 && <Scatter svgRef={refContext} xLabel="Largest healthy resting drift (|log2 FC|)" yLabel="Your difference (|log2 FC|)"
          points={context.withScore.map((o) => ({ x: o.drift!, y: Math.abs(molecules[o.i].score!), label: molecules[o.i].name }))} />}
      </div>
    </Card>),
    bars: (barItems.length > 0 && <Card title={ranked ? "Does MoTrPAC rank genes like your list?" : "Does the whole set shift?"}
      chips={[tissueLabel(f.tissue), f.contrast === ALL ? "All contrasts" : CATEGORY[f.contrast], cameraNote]}
      foot={<><button type="button" className="lr-ghost" onClick={() => csvDownload(`${layer.id}_comparisons.csv`, cols.map((c) => ({
        species: c.species, tissue: c.tissue_label, layer: c.layer, contrast: c.category, time: c.time_label, sex: c.sex, n_measured: c.n_measured,
        set_t: c.set_t, set_p: c.set_p, set_bh: c.set_bh, set_bonferroni: c.set_bonferroni, n_same: c.n_same, n_opposite: c.n_opposite,
        n_up: c.n_up, n_down: c.n_down, rho: c.rho, rho_p: c.rho_p, rho_bh: c.rho_bh, rho_bonferroni: c.rho_bonferroni })))}>CSV of all comparisons</button>
        <PngButton svgRef={refBars} name={`${layer.id}_set.png`} /></>}>
      <Bars svgRef={refBars} unit={ranked ? "Spearman ρ" : "cameraPR t"} labels={setLabels} items={barItems} />
    </Card>),
  };
  // Long lists lead with the set-level views; only the first rows fit in the heat map.
  const order = isMetab && human ? ["context", "molecules", "time", "overview", "bars"]
    : large ? ["overview", "bars", "molecules", "time"] : ["molecules", "time", "overview", "bars"];

  return <>
    <div className="xp-filters">
      <div className="lr-switch xp-species" role="group" aria-label="Species">
        {["human", "rat"].map((s) => <button type="button" key={s} aria-pressed={f.species === s} disabled={!cols.some((c) => c.species === s)}
          title={cols.some((c) => c.species === s) ? undefined : "Not measured in this layer"} onClick={() => pick("species", s)}>
          {s === "human" ? "Human · one exercise bout" : "Rat · 1–8 weeks training"}</button>)}
      </div>
      <Pick label="Tissue" value={f.tissue} onChange={(v) => pick("tissue", v)} options={tissues.map((t) => ({ value: t, label: tissueLabel(t) }))} />
      <Pick label="Contrast" value={f.contrast} onChange={(v) => pick("contrast", v)}
        options={[{ value: ALL, label: "All" }, ...CATEGORY_ORDER.filter((c) => c !== "TRAIN-SED" || !human).map((c) => ({
          value: c, label: CATEGORY[c], disabled: !categories.includes(c), note: categories.includes(c) ? undefined : human ? "Not measured in this tissue" : "Rats: training only" }))]} />
      {!human && <Pick label="Sex" value={f.sex} onChange={(v) => pick("sex", v)} options={[{ value: ALL, label: "Both" }, { value: "female", label: "Female" }, { value: "male", label: "Male" }]} />}
      <Pick label="Time" value={f.time} onChange={(v) => pick("time", v)} options={[{ value: ALL, label: "All" }, ...timesOf(inTissue).map((t) => ({ value: t, label: timeLabel(t) }))]} />
    </div>

    {order.map((key) => cards[key] ? <div key={key}>{cards[key]}</div> : null)}
  </>;
}

// ---- pathways ------------------------------------------------------------------------------------------

function PathwayView({ layer, rule }: { layer: ExplorerLayer; rule: Significance }) {
  const cols = layer.columns ?? [], values = layer.values ?? [], pathways = layer.molecules ?? [];
  const tissues = uniq(cols.map((c) => c.tissue));
  const [f, setF] = useState<Filters>({ species: "human", tissue: tissues.includes("blood") ? "blood" : tissues[0], layer: "RNA", contrast: "VS-CON", sex: ALL, time: ALL });
  const ref = useRef<SVGSVGElement>(null);
  const layers = uniq(cols.filter((c) => c.tissue === f.tissue).map((c) => c.layer));
  const selected = cols.filter((c) => c.tissue === f.tissue && c.layer === f.layer && inContrast(c, f.contrast));
  const index = new Map(cols.map((c, j) => [c, j]));
  const valueAt = new Map(values.map((v) => [`${v[0]}|${v[1]}`, v]));
  const groups = uniq(selected.map((c) => c.category)).sort((a, b) => rank(CATEGORY_ORDER, a) - rank(CATEGORY_ORDER, b)).map((cat) => ({
    title: CATEGORY[cat], cols: selected.filter((c) => c.category === cat).sort((a, b) => a.time_rank - b.time_rank).map((c) => ({ key: String(index.get(c)), label: shortTime(c.time_label) })) }));
  const cells = selected.flatMap((c) => pathways.map((_, i) => valueAt.get(`${i}|${index.get(c)}`)).filter(Boolean) as ExplorerValue[]);
  const up = cells.filter((v) => (v[2] ?? 0) > 0 && passes(valueP(v, rule))).length;
  const down = cells.filter((v) => (v[2] ?? 0) < 0 && passes(valueP(v, rule))).length;
  const positive = cells.filter((v) => (v[2] ?? 0) > 0).length;
  const label = (i: number) => `${pathways[i].name}${pathways[i].dir > 0 ? "  ↑" : pathways[i].dir < 0 ? "  ↓" : ""}`;
  const tissueLabel = (t: string) => cols.find((c) => c.tissue === t)?.tissue_label ?? t;
  return <>
    <div className="xp-filters">
      <Pick label="Tissue" value={f.tissue} onChange={(v) => setF({ ...f, tissue: v, layer: cols.some((c) => c.tissue === v && c.layer === f.layer) ? f.layer : cols.find((c) => c.tissue === v)!.layer })}
        options={tissues.map((t) => ({ value: t, label: tissueLabel(t) }))} />
      <Pick label="Layer" value={f.layer} onChange={(v) => setF({ ...f, layer: v })} options={layers.map((l) => ({ value: l, label: l }))} />
      <Pick label="Contrast" value={f.contrast} onChange={(v) => setF({ ...f, contrast: v })}
        options={[{ value: "VS-CON", label: CATEGORY["VS-CON"] }, { value: ALL, label: "All" }, ...CATEGORY_ORDER.filter((c) => c !== "TRAIN-SED").map((c) => ({ value: c, label: CATEGORY[c] }))]} />
    </div>
    <Card title={`Your pathways in healthy ${tissueLabel(f.tissue).toLowerCase()} after exercise`}
      chips={[`${f.layer} · MoTrPAC pathway tests (CAMERA)`, `${pathways.length} pathways matched`, "Human · one exercise bout", "↑ ↓ = direction in your list"]}
      stat={`${up}/${cells.length}`} caption={`pathway × time points rise significantly · ${down} fall · ${positive} point upward`}
      foot={<><span className="lr-hint">{layer.molecule_rule}</span>
        <button type="button" className="lr-ghost" onClick={() => csvDownload("pathways.csv", values.map((v) => {
          const c = cols[v[1]];
          return { pathway: pathways[v[0]].name, motrpac_set: pathways[v[0]].set, tissue: c.tissue_label, layer: c.layer, contrast: c.category, time: c.time_label, camera_z: v[2], raw_p: v[3], bh_q: v[4], bonferroni: v[5] };
        }))}>CSV</button><PngButton svgRef={ref} name="pathways_dotplot.png" /></>}>
      {groups.length ? <DotPlot svgRef={ref} rows={pathways.map((_, i) => label(i))} groups={groups} labels={["lower after exercise", "higher"]}
        cell={(row, key) => { const i = pathways.findIndex((_, k) => label(k) === row); const v = valueAt.get(`${i}|${key}`);
          return v ? { value: v[2], q: v[4], pass: passes(valueP(v, rule)), title: `${pathways[i].name} · ${cols[Number(key)].time_label}\nCAMERA z ${fmt(v[2])} · raw P ${fmt(v[3])} · BH q ${fmt(v[4])} · Bonferroni ${fmt(v[5])}` } : null; }} />
        : <div className="lr-empty">No MoTrPAC pathway test for this tissue, layer and contrast.</div>}
    </Card>
  </>;
}

// ---- page -------------------------------------------------------------------------------------------

export default function MotrpacExplorer() {
  const [examples, setExamples] = useState<ExplorerExample[]>([]);
  const [drafts, setDrafts] = useState<Record<ExplorerKind, Draft>>({ genes: EMPTY, metabolites: EMPTY, pathways: EMPTY });
  const [result, setResult] = useState<Load<ExplorerResponse>>({ kind: "idle" });
  const [rule, setRule] = useState<Significance>("bh");
  const [tab, setTab] = useState(0);
  const resultsRef = useRef<HTMLDivElement>(null);

  // Examples and their saved results ship with the website, so they work without the API.
  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}examples/index.json`).then((r) => r.ok ? r.json() : Promise.reject(r.status))
      .then((r: { examples: ExplorerExample[] }) => setExamples(r.examples))
      .catch(() => api.explorerExamples().then((r) => setExamples(r.examples)).catch(() => undefined));
  }, []);
  const show = (data: ExplorerResponse) => {
    setResult({ kind: "ready", data });
    setTab(Math.max(0, data.layers.findIndex((l) => l.available)));
    window.setTimeout(() => resultsRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }), 80);
  };
  const run = () => {
    const ready = KINDS.filter((k) => drafts[k.id].example || drafts[k.id].text.trim());
    const lists = ready.map((k) => drafts[k.id].example
      ? { kind: k.id, example: drafts[k.id].example!.id } : { kind: k.id, text: drafts[k.id].text, name: drafts[k.id].name ?? `Your ${k.label.toLowerCase()}` });
    setResult({ kind: "loading" });
    const live = () => api.explorerAnalyse({ lists }).then(show).catch((e) => setResult({ kind: "error", message: errorText(e) }));
    const example = ready.length === 1 ? drafts[ready[0].id].example : null;
    if (!example) { live(); return; }
    // An example loads its saved result; only if that file is missing does it fall back to the live analysis.
    fetch(`${import.meta.env.BASE_URL}examples/${example.id}.json`).then((r) => r.ok ? r.json() : Promise.reject(r.status))
      .then((data: ExplorerResponse) => show(data)).catch(live);
  };
  const data = result.kind === "ready" ? result.data : null;
  const layer = data?.layers[tab];
  const input = layer && layer.input !== null ? data!.inputs[layer.input] : null;

  return <section className="view lr">
    <header className="lr-hero">
      <div><span className="lr-kicker">MoTrPAC explorer</span><h2>Any molecule set, every exercise result</h2></div>
      <RuleSwitch rule={rule} onChange={setRule} />
    </header>
    <InputPanel examples={examples} drafts={drafts} setDrafts={setDrafts} onRun={run} running={result.kind === "loading"} />
    {result.kind === "error" && <div className="lr-empty">{result.message}</div>}
    {data && <div ref={resultsRef}>
      <div className="xp-inputs">{data.inputs.map((i) => <div key={i.id} className="xp-input-summary">
        <strong>{i.name}<em className={"xp-source" + (data.saved ? " is-saved" : "")}>{data.saved ? `Saved result · ${data.saved.generated_at}` : "Live analysis"}</em></strong>
        <span>{i.n_mapped} of {i.n_rows} matched{i.n_names && i.n_names !== i.n_mapped ? ` (${i.n_names} MoTrPAC names)` : ""}{i.directed ? " · with direction" : ""}{i.ranked ? " · ranked" : ""}</span>
        {i.unmapped.length > 0 && <details><summary>{i.unmapped.length}{i.unmapped.length === 50 ? "+" : ""} not found in MoTrPAC</summary><p>{i.unmapped.join(" · ")}</p></details>}
      </div>)}</div>
      <div className="lr-bar">
        <nav className="lr-tabs" aria-label="Omic layers">
          {data.layers.map((l, i) => <button type="button" key={`${l.id}-${i}`} disabled={!l.available} title={l.available ? undefined : l.reason}
            className={"lr-tab" + (tab === i ? " is-active" : "") + (!l.available ? " is-off" : "")} onClick={() => l.available && setTab(i)}>
            {l.label}{l.available ? <span>{l.n_matched > 999 ? `${Math.round(l.n_matched / 1000)}k` : l.n_matched}</span> : <em>{l.id === "epigenomics" ? "not loaded" : "no match"}</em>}
          </button>)}
        </nav>
      </div>
      <div className="lr-cards">
        {layer?.available ? (layer.id === "pathways" ? <PathwayView key={`p-${tab}`} layer={layer} rule={rule} />
          : <LayerView key={`${layer.id}-${tab}`} layer={layer} directed={!!input?.directed} ranked={!!input?.ranked} rule={rule} hint={input?.name} />)
          : <div className="lr-empty">{layer?.reason}</div>}
      </div>
      {input?.source && <p className="lr-footnote">Example source: {input.source}.</p>}
      <p className="lr-footnote">MoTrPAC values are published summary statistics (human: MotrpacHumanPreSuspensionAnalysis 2.0.8; rat: MotrpacRatTraining6moData 2.0.0). MoTrPAC participants and rats are healthy; comparisons with other studies describe association, not treatment effects.</p>
    </div>}
  </section>;
}
