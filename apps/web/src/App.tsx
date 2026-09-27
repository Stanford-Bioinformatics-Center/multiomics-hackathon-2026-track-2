import { createElement, useEffect, useState, type ReactNode } from "react";
import TechnicalFlowDiagram from "./components/TechnicalFlowDiagram";
import MotrpacExplorer from "./components/MotrpacExplorer";
import MetabolomicsCaseStudy from "./components/MetabolomicsCaseStudy";
import GeneralizedQuery from "./components/GeneralizedQuery";

type View = "architecture" | "technical" | "workflow" | "dashboard" | "live" | "metabolomics" | "generalized" | "slide";

const views: { id: View; label: string; eyebrow: string }[] = [
  { id: "architecture", label: "System architecture", eyebrow: "01" },
  { id: "technical", label: "Technical flow", eyebrow: "02" },
  { id: "workflow", label: "Researcher workflow", eyebrow: "03" },
  { id: "dashboard", label: "Reading guide", eyebrow: "04" },
  { id: "live", label: "Explorer", eyebrow: "05" },
  { id: "metabolomics", label: "Metabolomics (ST000763)", eyebrow: "06" },
  { id: "generalized", label: "Generalized query", eyebrow: "07" },
  { id: "slide", label: "Judge slide", eyebrow: "08" },
];

function viewFromUrl(): View {
  if (typeof window === "undefined") return "architecture";
  const section = new URLSearchParams(window.location.search).get("section");
  return views.find((item) => item.id === section)?.id ?? "architecture";
}

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
      "Cross-omic agreement labels",
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
  ["Label", "Describe disease and exercise directions, then RNA/protein agreement.", "classify"],
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
            <span className="pill pill--teal">Opposite direction from the PAH change</span>
            <span className="pill pill--coral">Same direction as the PAH change</span>
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
          <Text className="eyebrow">Cross-layer question</Text>
          <Text as="h3">Do RNA and protein move together in the same tissue and time window?</Text>
          <Text as="p">The PAH and MoTrPAC cohorts are separate. Agreement or disagreement between MoTrPAC layers describes healthy exercise physiology; it is a hypothesis for future patient studies.</Text>
        </div>
      </div>
    </section>
  );
}


function Dashboard({ onOpenLive }: { onOpenLive: () => void }) {
  return (
    <section className="view">
      <SectionIntro
        number="04"
        kicker="Reading guide"
        title="Read each contrast first"
        copy="The next section shows measured results. This guide explains how to read the panels without mixing the disease and exercise cohorts."
      />
      <div className="reading-comparison">
        <article>
          <span>REST · DISEASE STUDY</span>
          <h3>PAH compared with its own controls</h3>
          <p>The disease direction comes from a separate PAH study. Its numerator, denominator, tissue, and layer remain visible with each result.</p>
        </article>
        <article>
          <span>EXERCISE · MOTRPAC</span>
          <h3>Healthy physiology</h3>
          <p>Human data follow one acute bout in sedentary adults, compared with a resting control group. Rat data follow endurance training over weeks.</p>
        </article>
      </div>
      <div className="reading-note">
        <strong>Interpretation</strong>
        <p>Compare directions only after matching species, tissue, layer, time, and contrast. A cross-cohort match or difference generates a question for a future patient study; it does not test treatment in PAH.</p>
        <Button className="run-comparison" onClick={onOpenLive}>Open live results →</Button>
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
    ["05", "Two-axis description", "Disease direction + omic agreement"],
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
          <div><span className="result-mark">✓</span><strong>Compare</strong><small>RNA + protein directions within matched tissue and time</small></div>
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
  const [view, setView] = useState<View>(viewFromUrl);

  useEffect(() => {
    const syncFromUrl = () => setView(viewFromUrl());
    window.addEventListener("popstate", syncFromUrl);
    return () => window.removeEventListener("popstate", syncFromUrl);
  }, []);

  function navigate(next: View) {
    if (next === view) return;
    const url = new URL(window.location.href);
    url.searchParams.set("section", next);
    window.history.pushState(null, "", url);
    setView(next);
  }

  return (
    <main className="app-shell">
      <header className="app-header">
        <div className="app-identity">
          <div className="app-logo"><span>ES</span></div>
          <div><Text as="h1">Exercise Signature Explorer</Text><Text as="p">MoTrPAC exercise results for any molecule set</Text></div>
        </div>
        <div className="science-disclaimer"><span>RESEARCH USE</span> Directional evidence for follow-up — not a treatment claim</div>
      </header>
      <nav className="view-nav" aria-label="Deliverables">
        {views.map(item => (
          <Button className={`view-nav__item ${view === item.id ? "view-nav__item--active" : ""}`} onClick={() => navigate(item.id)} aria-current={view === item.id ? "page" : undefined} key={item.id}>
            <span>{item.eyebrow}</span>{item.label}
          </Button>
        ))}
      </nav>
      {view === "architecture" && <Architecture />}
      {view === "technical" && <TechnicalFlow />}
      {view === "workflow" && <Workflow />}
      {view === "dashboard" && <Dashboard onOpenLive={() => navigate("live")} />}
      {view === "live" && <MotrpacExplorer />}
      {view === "metabolomics" && <MetabolomicsCaseStudy />}
      {view === "generalized" && <GeneralizedQuery />}
      {view === "slide" && <JudgeSlide />}
      <footer className="app-footer"><span>System handoff · v0.1</span><span>Cross-cohort comparison · transparent assumptions · reproducible export</span></footer>
    </main>
  );
}
