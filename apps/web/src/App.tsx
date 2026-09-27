import { createElement, useState, type ReactNode } from "react";
import QueryBuilder, { initialAnalysisQuery } from "./components/QueryBuilder";
import TechnicalFlowDiagram from "./components/TechnicalFlowDiagram";
import { getStudyContext, getVisualizationMode, type AnalysisQuery } from "./domain/analysis";

type View = "architecture" | "technical" | "workflow" | "dashboard" | "slide";
type DashboardTab = "Disease vs Exercise" | "RNA vs Protein" | "Heatmap" | "Evidence Table" | "Methods & Limitations";

const views: { id: View; label: string; eyebrow: string }[] = [
  { id: "architecture", label: "System architecture", eyebrow: "01" },
  { id: "technical", label: "Technical flow", eyebrow: "02" },
  { id: "workflow", label: "Researcher workflow", eyebrow: "03" },
  { id: "dashboard", label: "Results dashboard", eyebrow: "04" },
  { id: "slide", label: "Judge slide", eyebrow: "05" },
];

const dashboardTabs: DashboardTab[] = [
  "Disease vs Exercise",
  "RNA vs Protein",
  "Heatmap",
  "Evidence Table",
  "Methods & Limitations",
];

function Text({
  as = "div",
  className = "",
  children,
  ...props
}: {
  as?: "div" | "span" | "p" | "h1" | "h2" | "h3";
  className?: string;
  children: ReactNode;
  [key: string]: unknown;
}) {
  return createElement(as, { className, ...props }, children);
}

function Button({
  className = "",
  children,
  ...props
}: {
  className?: string;
  children: ReactNode;
  [key: string]: unknown;
}) {
  return createElement("button", { className, type: "button", ...props }, children);
}

function Arrow({ vertical = false }: { vertical?: boolean }) {
  return <span className={vertical ? "arrow arrow--vertical" : "arrow"} aria-hidden="true" />;
}

function Mark({ type = "solid" }: { type?: "solid" | "planned" | "assumption" }) {
  return (
    <span className={`status-mark status-mark--${type}`}>
      {type === "planned" ? "PLANNED" : type === "assumption" ? "ASSUMPTION" : "PAH DEMO"}
    </span>
  );
}

function MiniNode({
  children,
  planned = false,
  badge,
}: {
  children: ReactNode;
  planned?: boolean;
  badge?: string;
}) {
  return (
    <div className={`mini-node ${planned ? "mini-node--planned" : ""}`}>
      <span>{children}</span>
      {badge && <span className="mini-node__badge">{badge}</span>}
    </div>
  );
}

function SectionIntro({
  number,
  kicker,
  title,
  copy,
}: {
  number: string;
  kicker: string;
  title: string;
  copy: string;
}) {
  return (
    <div className="section-intro">
      <div className="section-intro__number">{number}</div>
      <div>
        <Text className="eyebrow">{kicker}</Text>
        <Text as="h2" className="section-title">{title}</Text>
        <Text as="p" className="section-copy">{copy}</Text>
      </div>
    </div>
  );
}

const architectureLayers = [
  {
    number: "01",
    title: "User + frontend",
    tone: "navy",
    nodes: [
      "Researcher",
      "Curated PAH example",
      "Structured CSV upload",
      "Query builder",
      "Context filters",
      "Results dashboard",
      "Export",
    ],
  },
  {
    number: "02",
    title: "Application API",
    tone: "navy",
    nodes: [
      "Signature validation",
      "Dataset availability resolver",
      "Query orchestration",
      "Results API",
      "Provenance API",
      "Export service",
    ],
  },
  {
    number: "03",
    title: "Scientific processing",
    tone: "purple",
    nodes: [
      "Identifier normalization",
      "Alias resolution",
      "Human → rat ortholog mapping",
      "Omic-specific adapters",
      "Disease ↔ exercise comparison",
      "Cross-omic classifier",
      "Pathway / context annotation",
      "Evidence + limitation labels",
    ],
  },
  {
    number: "04",
    title: "Data products",
    tone: "orange",
    nodes: [
      "Directed signature registry",
      "MoTrPAC transcriptomics",
      "MoTrPAC proteomics",
      "Identifier + ortholog tables",
      "Pathway resources",
      "Versioned manifests",
      "Study-context catalog",
      "Layer capability registry",
    ],
  },
  {
    number: "05",
    title: "Outputs",
    tone: "teal",
    nodes: [
      "Disease ↔ exercise plot",
      "RNA ↔ protein plot",
      "Gene-by-layer heatmap",
      "Evidence summary cards",
      "Filterable results",
      "Mapping + provenance report",
      "CSV / JSON bundle",
    ],
  },
];

function Architecture() {
  return (
    <section className="view">
      <SectionIntro
        number="01"
        kicker="Technical handoff"
        title="A deterministic evidence pipeline"
        copy="Five inspectable layers carry a directed disease signature into a cross-cohort, multi-omic comparison. Every transform remains visible."
      />

      <div className="architecture-shell">
        <div className="architecture-main">
          {architectureLayers.map((layer, index) => (
            <div className="architecture-row-wrap" key={layer.title}>
              <div className={`architecture-row architecture-row--${layer.tone}`}>
                <div className="architecture-label">
                  <span className="architecture-label__number">{layer.number}</span>
                  <Text as="h3">{layer.title}</Text>
                </div>
                <div className="architecture-nodes">
                  {layer.nodes.map((node) => (
                    <MiniNode key={node} badge={node === "Context filters" ? "source species · target species · derived design · tissue · sex · time point · omics" : undefined}>
                      {node}
                    </MiniNode>
                  ))}
                  {layer.number === "01" && <MiniNode planned>Natural-language query</MiniNode>}
                  {layer.number === "04" && <MiniNode planned>Genomics adapter</MiniNode>}
                  {layer.number === "04" && <MiniNode planned>Epigenomics adapter</MiniNode>}
                  {layer.number === "04" && <MiniNode planned>Metabolomics adapter</MiniNode>}
                  {layer.number === "05" && <MiniNode planned>Predictive modeling</MiniNode>}
                </div>
              </div>
              {index < architectureLayers.length - 1 && <Arrow vertical />}
            </div>
          ))}
        </div>

        <aside className="llm-rail">
          <div className="llm-card">
            <Mark type="planned" />
            <Text as="h3">Optional LLM assistant</Text>
            <Text as="p">Convenience layer, isolated from the scientific calculation path.</Text>
            <div className="llm-link"><span>01</span> Natural language → structured query</div>
            <div className="llm-link"><span>02</span> Plain-language result summaries</div>
          </div>
          <div className="legend-card">
            <Text className="eyebrow">System key</Text>
            <div className="legend-row"><span className="legend-swatch legend-swatch--solid" /> Available in PAH demo</div>
            <div className="legend-row"><span className="legend-swatch legend-swatch--dashed" /> Planned capability</div>
            <div className="legend-row"><span className="legend-swatch legend-swatch--warning" /> Assumption / caution</div>
          </div>
        </aside>
      </div>

      <div className="deterministic-note">
        <span className="deterministic-note__icon">✓</span>
        <Text><strong>Deterministic by design.</strong> All identifiers, contrasts, mappings, and calculations are deterministically validated.</Text>
      </div>
      <div className="assumption-strip">
        <Mark type="assumption" />
        <Text><strong>Human → rat mapping</strong> and <strong>RNA → protein comparison</strong> are explicit analytical assumptions, preserved in every provenance report.</Text>
      </div>
    </section>
  );
}

function TechnicalFlow() {
  return (
    <section className="view">
      <SectionIntro
        number="02"
        kicker="Deployment and data flow"
        title="Every interaction has a named contract"
        copy="The target species resolves a read-only study design from the catalog before availability, mapping, or scientific processing begins."
      />
      <TechnicalFlowDiagram />
      <div className="api-contract">
        <span>GET /api/v1/catalog</span>
        <span>GET /api/v1/catalog/species/:species</span>
        <span>POST /api/v1/signatures/validate</span>
        <span>POST /api/v1/availability</span>
        <span>POST /api/v1/comparisons</span>
        <span>GET /api/v1/comparisons/:id/features</span>
        <span>GET /api/v1/comparisons/:id/provenance</span>
        <span>GET /api/v1/comparisons/:id/export</span>
      </div>
      <div className="flow-contract-note">
        <span>VALIDATION CONTRACT</span>
        <strong>All identifiers, contexts, mappings, and calculations are deterministic.</strong>
        <p>The optional LLM may propose a query, explain an empty state, or summarize validated output. It cannot alter species, omics, identifiers, effects, significance, or classifications.</p>
      </div>
    </section>
  );
}

const workflowSteps = [
  ["Choose", "Select the curated PAH example or upload a directed signature.", "signature"],
  ["Review", "Inspect recognized identifiers, directions, effects, and source.", "input"],
  ["Resolve", "Resolve aliases and human-to-rat species mappings.", "mapping"],
  ["Confirm", "Approve, edit, or exclude every normalization decision.", "checkpoint"],
  ["Resolve context", "Choose the target MoTrPAC species; study design is derived from the catalog.", "context"],
  ["Compare", "Choose compatible omics, tissue, sex, and time point, then run precomputed results.", "compute"],
  ["Classify", "Label disease opposition and RNA/protein agreement.", "classify"],
  ["Explore", "Review plots, pathways, provenance, and limitations.", "evidence"],
  ["Export", "Download a reproducible CSV/JSON result bundle.", "export"],
];

function Workflow() {
  return (
    <section className="view">
      <SectionIntro
        number="02"
        kicker="Researcher workflow"
        title="No feature disappears in the pipeline"
        copy="A deliberate confirmation checkpoint separates computational normalization from researcher-approved analysis."
      />
      <div className="workflow-context">
        <div>
          <Text className="eyebrow">Current comparison contract</Text>
          <Text as="h3">PAH skeletal muscle <span>vs healthy control</span></Text>
        </div>
        <Arrow />
        <div>
          <Text className="eyebrow">Exercise response</Text>
          <Text as="h3">Chronically trained rat <span>vs sex-matched sedentary rat</span></Text>
        </div>
      </div>
      <div className="workflow-track">
        {workflowSteps.map(([title, copy, type], index) => (
          <div className={`workflow-step ${type === "checkpoint" ? "workflow-step--checkpoint" : ""}`} key={title}>
            <div className="workflow-step__top">
              <span className="workflow-step__number">{String(index + 1).padStart(2, "0")}</span>
              <span className="workflow-step__line" />
            </div>
            <div className="workflow-step__card">
              {type === "checkpoint" && <span className="checkpoint-label">REQUIRED CONFIRMATION</span>}
              <Text as="h3">{title}</Text>
              <Text as="p">{copy}</Text>
            </div>
          </div>
        ))}
      </div>
      <div className="workflow-safety">
        <div className="workflow-safety__mark">!</div>
        <div>
          <Text as="h3">Explicit mapping contract</Text>
          <Text as="p">The application must never silently drop or remap a feature. Unmapped and ambiguous records remain visible, downloadable, and attributable to the source row.</Text>
        </div>
      </div>
      <div className="workflow-availability">
        <div><Text className="eyebrow">Catalog-derived context</Text><Text as="h3">The researcher selects species—not study design</Text></div>
        <div className="availability-example"><span>Target species</span><strong>Human</strong><small>Acute exercise study</small></div>
        <div className="availability-example"><span>Target species</span><strong>Rat</strong><small>Chronic exercise-training study</small></div>
        <div className="availability-example availability-example--empty"><span>Context mismatch</span><strong>Scientific empty state</strong><small>Alternatives remain visible</small></div>
        <div className="availability-example availability-example--empty"><span>Planned layer</span><strong>Unsupported omic</strong><small>Run supported layers instead</small></div>
      </div>
      <div className="classification-grid">
        <div className="classification-card">
          <Text className="eyebrow">Independent classification A</Text>
          <Text as="h3">Disease relationship</Text>
          <div className="classification-list">
            <span className="pill pill--teal">Opposes disease direction</span>
            <span className="pill pill--coral">Matches disease direction</span>
            <span className="pill">Near-zero / uncertain</span>
            <span className="pill">Missing / unmapped</span>
          </div>
        </div>
        <div className="classification-card">
          <Text className="eyebrow">Independent classification B</Text>
          <Text as="h3">Cross-omic relationship — MVP</Text>
          <div className="classification-list">
            <span className="pill pill--teal">Concordant up / concordant down</span>
            <span className="pill pill--coral">Higher measured RNA · lower measured protein</span>
            <span className="pill pill--coral">Lower measured RNA · higher measured protein</span>
            <span className="pill">RNA only</span>
            <span className="pill">Protein only</span>
            <span className="pill">Near-zero</span>
            <span className="pill">Insufficient data</span>
          </div>
        </div>
        <div className="strongest-state">
          <Text className="eyebrow">Strongest exploratory state</Text>
          <Text as="h3">RNA and protein agree — and both oppose the disease direction.</Text>
          <Text as="p">A prioritized signal for hypothesis generation and follow-up, not evidence that exercise treats PAH. RNA and protein remain independent measurements; neither is treated as a proxy for the other.</Text>
        </div>
      </div>
    </section>
  );
}

const genes = [
  { name: "PPARGC1A", x: 63, y: 29, type: "oppose" },
  { name: "SOD2", x: 72, y: 24, type: "oppose" },
  { name: "NDUFS1", x: 67, y: 37, type: "oppose" },
  { name: "COL1A1", x: 29, y: 70, type: "oppose" },
  { name: "CXCL12", x: 31, y: 61, type: "oppose" },
  { name: "STAT3", x: 69, y: 68, type: "align" },
  { name: "IL6R", x: 76, y: 76, type: "align" },
  { name: "MYH7", x: 26, y: 26, type: "align" },
  { name: "VEGFA", x: 57, y: 47, type: "neutral" },
];

function QuadrantPlot({ omic = false }: { omic?: boolean }) {
  return (
    <div className="plot-shell">
      <div className="plot-y-label">{omic ? "Protein log₂FC" : "Exercise log₂FC"}</div>
      <div className="plot">
        <div className="quadrant quadrant--tl"><span>{omic ? "Discordant" : "Opposes"}</span></div>
        <div className="quadrant quadrant--tr"><span>{omic ? "Concordant +" : "Disease aligned"}</span></div>
        <div className="quadrant quadrant--bl"><span>{omic ? "Concordant −" : "Disease aligned"}</span></div>
        <div className="quadrant quadrant--br"><span>{omic ? "Discordant" : "Opposes"}</span></div>
        <span className="axis axis--x" />
        <span className="axis axis--y" />
        {genes.map((gene, index) => (
          <span
            className={`plot-point plot-point--${gene.type}`}
            style={{ left: `${gene.x}%`, top: `${gene.y}%` }}
            key={gene.name}
          >
            {(index < 5 || omic) && <span className="plot-label">{gene.name}</span>}
          </span>
        ))}
      </div>
      <div className="plot-x-label">{omic ? "RNA log₂FC" : "Disease effect"}</div>
    </div>
  );
}

function Heatmap() {
  const rows = ["PPARGC1A", "SOD2", "NDUFS1", "COL1A1", "CXCL12", "STAT3", "IL6R", "MYH7"];
  const cells = [
    ["down", "up", "up", "na"], ["down", "up", "up", "na"], ["down", "up", "up", "na"],
    ["up", "down", "down", "na"], ["up", "down", "down", "na"], ["up", "up", "up", "na"],
    ["up", "up", "down", "na"], ["down", "down", "down", "na"],
  ];
  return (
    <div className="heatmap">
      <div className="heatmap__head"><span />{["Disease", "RNA", "Protein", "Metabolite"].map((h) => <span key={h}>{h}</span>)}</div>
      {rows.map((row, i) => (
        <div className="heatmap__row" key={row}>
          <strong>{row}</strong>
          {cells[i].map((cell, j) => <span className={`heat-cell heat-cell--${cell}`} key={j}>{j !== 3 && (i + j) % 3 === 0 ? "•" : ""}</span>)}
        </div>
      ))}
      <div className="heatmap__legend"><span className="legend-gradient" /> signed effect <b>•</b> adjusted significant <i>planned layer</i></div>
    </div>
  );
}

function EvidenceTable() {
  const rows = [
    ["PPARGC1A", "PPARGC1A → Ppargc1a", "−1.24", "+0.88", "+0.61", "Opposes · concordant"],
    ["SOD2", "SOD2 → Sod2", "−0.76", "+0.54", "+0.42", "Opposes · concordant"],
    ["COL1A1", "COL1A1 → Col1a1", "+1.45", "−0.63", "−0.38", "Opposes · concordant"],
    ["STAT3", "STAT3 → Stat3", "+0.71", "+0.39", "+0.26", "Matches · concordant"],
    ["VEGFA", "VEGFA → Vegfa", "+0.32", "+0.05", "—", "Near-zero · RNA only"],
  ];
  return (
    <div className="evidence-table">
      <div className="table-tools"><div className="search-box">Search gene, mapping, source…</div><Button className="tiny-button">Download visible rows</Button></div>
      <div className="table-row table-row--head">{["Feature", "Mapping", "Disease", "RNA", "Protein", "Classification"].map(h => <span key={h}>{h}</span>)}</div>
      {rows.map(row => <div className="table-row" key={row[0]}>{row.map((cell, i) => <span key={i} className={i === 5 ? (cell.startsWith("Opposes") ? "good-text" : "warn-text") : ""}>{cell}</span>)}</div>)}
      <div className="table-provenance">＋ Expandable provenance rows include source accession, numerator, denominator, q-values, ortholog method, mapping confidence, and pipeline version.</div>
    </div>
  );
}

function Methods() {
  return (
    <>
      <div className="methods-grid">
        <div className="method-block"><Text className="eyebrow">Source studies</Text><Text as="h3">Cross-cohort by design</Text><Text as="p">External human PAH skeletal-muscle signature compared with a catalog-resolved MoTrPAC context.</Text></div>
        <div className="method-block"><Text className="eyebrow">Versions</Text><Text as="h3">Reproducible manifests</Text><Text as="p">Every export records dataset release, package versions, pipeline commit, parameters, timestamp, and content checksums.</Text></div>
        <div className="method-block"><Text className="eyebrow">Mappings</Text><Text as="h3">Assumptions stay visible</Text><Text as="p">Human-to-rat orthology and transcript-to-protein comparison are reported with method, confidence, ambiguity, and exclusions.</Text></div>
        <div className="method-block method-block--warning"><Text className="eyebrow">What is not supported</Text><Text as="h3">No treatment or causal claim</Text><Text as="p">Directional opposition is exploratory evidence for follow-up. It does not show that exercise treats PAH or establish clinical benefit.</Text></div>
      </div>
      <div className="scientific-state-gallery">
        <div className="state-gallery-card"><span>SCIENTIFIC EMPTY STATE</span><strong>No matching tissue or time point</strong><p>The resolved study exists, but no result matches the selected catalog filters. Show compatible alternatives.</p></div>
        <div className="state-gallery-card"><span>SCIENTIFIC EMPTY STATE</span><strong>No mapped features</strong><p>Show submitted identifiers, namespaces, alias attempts, mapping failures, and ortholog status.</p></div>
        <div className="state-gallery-card state-gallery-card--weak"><span>MEASUREMENTS RETAINED</span><strong>No adjusted-significant support</strong><p>Effect sizes, directions, and adjusted p-values remain visible rather than being treated as missing.</p></div>
        <div className="state-gallery-card state-gallery-card--error"><span>TECHNICAL ERROR ONLY</span><strong>Processing failure</strong><p>Reserved for malformed files, unavailable services, corrupted artifacts, or processing exceptions.</p></div>
      </div>
    </>
  );
}

function PairwiseMatrix() {
  return (
    <div className="pairwise-matrix">
      <div className="pairwise-matrix__head"><span /><strong>RNA</strong><strong>Protein</strong><strong>Metabolite</strong></div>
      {[
        ["RNA", "—", "0.68", "planned"],
        ["Protein", "0.68", "—", "planned"],
        ["Metabolite", "planned", "planned", "—"],
      ].map((row) => (
        <div className="pairwise-matrix__row" key={row[0]}>
          <strong>{row[0]}</strong>
          {row.slice(1).map((cell, index) => <span className={cell === "planned" ? "matrix-planned" : ""} key={index}>{cell}</span>)}
        </div>
      ))}
      <small>Pairwise concordance matrix · select a supported pair for detailed inspection</small>
    </div>
  );
}

function DashboardContent({ tab, query }: { tab: DashboardTab; query: AnalysisQuery }) {
  const mode = getVisualizationMode(query.selectedOmics);
  if (tab === "Disease vs Exercise") return <QuadrantPlot />;
  if (tab === "RNA vs Protein") {
    if (mode === "single_layer") {
      return (
        <div className="visual-empty-state">
          <span>ONE LAYER SELECTED</span>
          <Text as="h3">Disease-versus-exercise analysis remains available.</Text>
          <Text as="p">Cross-omic concordance or discordance requires at least two compatible layers. No missing biological mechanism is inferred.</Text>
        </div>
      );
    }
    if (mode === "matrix") {
      return <div className="multi-view"><Heatmap /><PairwiseMatrix /></div>;
    }
    return <QuadrantPlot omic />;
  }
  if (tab === "Heatmap") return <Heatmap />;
  if (tab === "Evidence Table") return <EvidenceTable />;
  return <Methods />;
}

function Dashboard() {
  const [tab, setTab] = useState<DashboardTab>("Disease vs Exercise");
  const [query, setQuery] = useState<AnalysisQuery>(initialAnalysisQuery);
  const context = getStudyContext(query.targetMotrpacSpecies)!;
  const isRat = query.targetMotrpacSpecies === "rat";
  const cards = [
    ["128", "Input features", "100% retained"],
    ["112", "Successfully mapped", "87.5%"],
    ["47", "Opposite direction", "42% of mapped"],
    ["31", "Adjusted-significant", "cross-cohort support"],
    ["26", "RNA / protein agree", "of 39 dual-measured"],
    ["16", "Missing / ambiguous", "requires review"],
  ];
  return (
    <section className="view">
      <SectionIntro number="03" kicker="Desktop wireframe" title="Evidence first, claims constrained" copy="The dashboard keeps contrast definitions, assumptions, missingness, and provenance beside the scientific signal." />
      <div className="dashboard-frame">
        <div className="dashboard-topbar">
          <div className="dashboard-brand"><span className="brand-mark">ES</span><div><strong>Exercise Signature Explorer</strong><small>Cross-cohort molecular comparison</small></div></div>
          <div className="topbar-context"><span>Current signature</span><strong>PAH skeletal muscle v1.2</strong></div>
          <div className="topbar-context"><span>Dataset / version</span><strong>{context.datasetId} · {context.datasetVersion}</strong></div>
          <Button className="method-link">Methodology ↗</Button>
          <Button className="export-button">Export bundle ↓</Button>
        </div>
        <div className="contrast-banner">
          <div><span>Disease numerator</span><strong>Human PAH skeletal muscle</strong></div><div className="versus">VS</div><div><span>Disease denominator</span><strong>Healthy human control</strong></div>
          <div className="contrast-divider" />
          <div><span>Exercise numerator</span><strong>{isRat ? "Chronically trained rat" : "Human post-acute exercise"}</strong></div><div className="versus">VS</div><div><span>Exercise denominator</span><strong>{isRat ? "Sex-matched sedentary rat" : "Matched human baseline"}</strong></div>
        </div>
        <div className="dashboard-body">
          <aside className="filter-rail">
            <QueryBuilder onQueryChange={setQuery} />
          </aside>
          <div className="dashboard-main">
            <div className="summary-cards">
              {cards.map(([value, label, detail], i) => <div className={`summary-card ${i === 2 || i === 4 ? "summary-card--focus" : ""}`} key={label}><strong>{value}</strong><span>{label}</span><small>{detail}</small></div>)}
            </div>
            <div className="visual-card">
              <div className="dashboard-tabs">
                {dashboardTabs.map(item => <Button className={`dashboard-tab ${tab === item ? "dashboard-tab--active" : ""}`} onClick={() => setTab(item)} key={item}>{item}</Button>)}
              </div>
              <div className="visual-header">
                <div><Text className="eyebrow">{tab === "Methods & Limitations" ? "Interpretation guardrails" : `${context.species} · ${context.studyDesign} · ${query.tissue} · ${query.sex} · ${query.timepoint}`}</Text><Text as="h3">{tab}</Text></div>
                {tab !== "Methods & Limitations" && <div className="chart-legend"><span><i className="dot dot--teal" /> Opposes / concordant</span><span><i className="dot dot--coral" /> Aligned / discordant</span><span><i className="dot dot--gray" /> Uncertain</span></div>}
              </div>
              <DashboardContent tab={tab} query={query} />
              {(tab === "Disease vs Exercise" || (tab === "RNA vs Protein" && query.selectedOmics.length === 2)) && <div className="tooltip-demo"><strong>PPARGC1A</strong><span>Human → rat: Ppargc1a · 1:1 ortholog</span><span>Disease effect: −1.24 · q 0.003</span><span>Exercise RNA: +0.88 · q 0.012</span><span>Source: PAH-EXT-04 / MoTrPAC 1.0</span></div>}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

function JudgeSlide() {
  const steps = [
    ["01", "Directed disease signature", "Signed effects + source"],
    ["02", "Species + omic selection", "Study design derived from catalog"],
    ["03", "Availability + mapping validation", "No silent reinterpretation"],
    ["04", "Compatible MoTrPAC data", "Independent RNA + protein"],
    ["05", "Two-axis classification", "Opposition + concordance"],
    ["06", "Interactive evidence dashboard", "Plots + provenance + export"],
  ];
  return (
    <section className="view">
      <SectionIntro number="04" kicker="Presentation slide" title="One contract, extensible evidence" copy="A concise version of the system story for review panels and judges." />
      <div className="judge-slide">
        <div className="judge-slide__head">
          <div><span className="judge-kicker">EXERCISE SIGNATURE EXPLORER</span><Text as="h2">From a directed signature to<br />inspectable multi-omic evidence.</Text></div>
          <div className="judge-thesis">Cross-cohort comparison<br /><strong>for hypothesis generation</strong></div>
        </div>
        <div className="judge-flow">
          {steps.map(([number, title, detail], index) => (
            <div className="judge-step-wrap" key={number}>
              <div className={`judge-step judge-step--${index + 1}`}><span>{number}</span><Text as="h3">{title}</Text><Text as="p">{detail}</Text></div>
              {index < steps.length - 1 && <Arrow />}
            </div>
          ))}
        </div>
        <div className="judge-result">
          <div><span className="result-mark">✓</span><strong>Prioritize</strong><small>Concordant RNA + protein that oppose disease direction</small></div>
          <div><span className="result-mark">↗</span><strong>Follow up</strong><small>Generate hypotheses; do not infer treatment or causality</small></div>
          <div><span className="result-mark">≡</span><strong>Reproduce</strong><small>Export mappings, versions, assumptions, and calculations</small></div>
        </div>
        <div className="judge-availability">
          <strong>Current MVP</strong>
          <span>Human → acute context from catalog</span>
          <span>Rat → chronic context from catalog</span>
          <span>Transcriptomics + proteomics implemented</span>
          <span>Additional omics use the same planned adapter contract</span>
        </div>
        <div className="judge-caption">“PAH skeletal muscle is our worked example; the same contract can accept additional directed signatures and omic-specific adapters.”</div>
      </div>
    </section>
  );
}

export default function App() {
  const [view, setView] = useState<View>("architecture");
  return (
    <main className="app-shell">
      <header className="app-header">
        <div className="app-identity">
          <div className="app-logo"><span>ES</span></div>
          <div><Text as="h1">Exercise Signature Explorer</Text><Text as="p">System-design handoff · PAH worked example</Text></div>
        </div>
        <div className="science-disclaimer"><span>RESEARCH USE</span> Directional evidence for follow-up — not a treatment claim</div>
      </header>
      <nav className="view-nav" aria-label="Deliverables">
        {views.map(item => (
          <Button className={`view-nav__item ${view === item.id ? "view-nav__item--active" : ""}`} onClick={() => setView(item.id)} key={item.id}>
            <span>{item.eyebrow}</span>{item.label}
          </Button>
        ))}
      </nav>
      {view === "architecture" && <Architecture />}
      {view === "technical" && <TechnicalFlow />}
      {view === "workflow" && <Workflow />}
      {view === "dashboard" && <Dashboard />}
      {view === "slide" && <JudgeSlide />}
      <footer className="app-footer"><span>System handoff · v0.1</span><span>Cross-cohort comparison · transparent assumptions · reproducible export</span></footer>
    </main>
  );
}
