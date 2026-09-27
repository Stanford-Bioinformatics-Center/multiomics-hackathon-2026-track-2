import { createElement, useEffect, useState, type ReactNode } from "react";
import LiveDashboard from "./components/LiveDashboard";
import MetabolomicsCaseStudy from "./components/MetabolomicsCaseStudy";
import GeneralizedQuery from "./components/GeneralizedQuery";
import Discordance from "./components/Discordance";

// ---------------------------------------------------------------------------
// Navigation model
// ---------------------------------------------------------------------------
//
// The application presents exactly five live views. The mock design views that
// used to live here (architecture/01, technical/02, workflow/03, dashboard/04,
// slide/08) have been retired — their content was migrated into
// DESIGN_PROVENANCE.md with the original Figma-Make attribution preserved.
//
// `View` therefore only names the five live views. Retired identifiers are kept
// in `RETIRED_VIEWS` purely so the redirect resolver can recognise a stale or
// hand-typed identifier and map it to a live view; they are never rendered.

export type View = "discordance" | "live" | "generalized" | "metabolomics" | "methods";

/** Identifiers of the removed mock design views. Never rendered. */
export const RETIRED_VIEWS = ["architecture", "technical", "workflow", "dashboard", "slide"] as const;

/** The default live view shown on initial load and used as the ultimate fallback. */
export const DEFAULT_VIEW: View = "discordance";

/** sessionStorage key holding the last live view the user visited this session. */
export const LAST_VIEW_STORAGE_KEY = "esx.lastView";

export const NAV_VIEWS: { id: View; label: string; eyebrow: string }[] = [
  { id: "discordance", label: "Discordance", eyebrow: "01" },
  { id: "live", label: "Live results (API)", eyebrow: "02" },
  { id: "generalized", label: "Generalized query", eyebrow: "03" },
  { id: "metabolomics", label: "Metabolomics (ST000763)", eyebrow: "04" },
  { id: "methods", label: "About / Methods", eyebrow: "05" },
];

const LIVE_VIEW_IDS = new Set<string>(NAV_VIEWS.map((entry) => entry.id));

/** True for any identifier that names a live, navigable view. */
export function isLiveView(candidate: string | null | undefined): candidate is View {
  return typeof candidate === "string" && LIVE_VIEW_IDS.has(candidate);
}

/**
 * Resolve any requested view identifier to a live view that is safe to render.
 *
 * - A live identifier resolves to itself.
 * - A retired or otherwise unknown identifier resolves to `lastLiveView` when that
 *   is itself a live view, otherwise to the default Discordance view.
 *
 * A retired view is NEVER returned, so callers can render the result directly.
 */
export function resolveView(requested: string | null | undefined, lastLiveView?: string | null): View {
  if (isLiveView(requested)) {
    return requested;
  }
  if (isLiveView(lastLiveView)) {
    return lastLiveView;
  }
  return DEFAULT_VIEW;
}

/**
 * Read any externally-requested view (a `?view=`/`#hash` value) plus the persisted
 * last-used live view from sessionStorage, and resolve them to the initial live view.
 * Runs before first render so a retired identifier is mapped to a live view up front.
 */
export function initialView(): View {
  let requested: string | null = null;
  let lastLiveView: string | null = null;

  if (typeof window !== "undefined") {
    try {
      const params = new URLSearchParams(window.location.search);
      requested = params.get("view") ?? (window.location.hash.replace(/^#\/?/, "") || null);
    } catch {
      requested = null;
    }
    try {
      lastLiveView = window.sessionStorage?.getItem(LAST_VIEW_STORAGE_KEY) ?? null;
    } catch {
      lastLiveView = null;
    }
  }

  return resolveView(requested, lastLiveView);
}

// ---------------------------------------------------------------------------
// Small presentational helpers reused by the About / Methods page
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// About / Methods — the single combined methods page (replaces the mock views)
// ---------------------------------------------------------------------------

const METHODS_DOCS: { href: string; title: string; blurb: string }[] = [
  {
    href: "../../DESIGN_PROVENANCE.md",
    title: "DESIGN_PROVENANCE.md",
    blurb:
      "Migrated content of the retired mock design views (system architecture, technical flow, researcher workflow, results dashboard, judge slide), with the original Figma-Make design attribution preserved.",
  },
  {
    href: "../../WORKFLOWS.md",
    title: "WORKFLOWS.md",
    blurb:
      "Traced map of every analysis workflow — entry point, module, data read, transform, output fields, provenance/run_id, and the UI view that renders each output.",
  },
  {
    href: "../../METABOLOMICS.md",
    title: "METABOLOMICS.md",
    blurb: "Methods and provenance for the ST000763 plasma metabolomics case study.",
  },
];

function AboutMethods() {
  return (
    <section className="view">
      <div className="section-intro">
        <div className="section-intro__number">05</div>
        <div>
          <Text className="eyebrow">About / Methods</Text>
          <Text as="h2" className="section-title">Honest framing and where the details live</Text>
          <Text as="p" className="section-copy">
            Reported results are cross-cohort association findings — they are not a claim of PAH
            treatment efficacy. Cohorts (muscle protein, plasma metabolomics, blood) are analyzed
            separately and never merged. The PTM occupancy caveat applies: a PTM signal is not
            automatically a measure of modification occupancy, enzyme activity, or functional
            consequence. The statistical engines remain the single source of truth; the API and this
            client only serve, shape, and render — they compute no statistics.
          </Text>
        </div>
      </div>

      <div className="methods-grid">
        {METHODS_DOCS.map((doc) => (
          <a className="method-block" href={doc.href} key={doc.href}>
            <Text className="eyebrow">Documentation</Text>
            <Text as="h3">{doc.title}</Text>
            <Text as="p">{doc.blurb}</Text>
          </a>
        ))}
      </div>

      <div className="assumption-strip">
        <span className="status-mark status-mark--assumption">DESIGN ORIGIN</span>
        <Text>
          <strong>Figma-Make design provenance.</strong> The earlier mock design views were produced
          in Figma-Make. Their content and attribution now live in{" "}
          <a href="../../DESIGN_PROVENANCE.md">DESIGN_PROVENANCE.md</a>.
        </Text>
      </div>
    </section>
  );
}

// ---------------------------------------------------------------------------
// App shell
// ---------------------------------------------------------------------------

export default function App() {
  const [view, setView] = useState<View>(initialView);

  // Persist the last-used live view so a later direct navigation to a retired
  // identifier can fall back to it (see `resolveView`/`initialView`).
  useEffect(() => {
    if (typeof window === "undefined") return;
    try {
      window.sessionStorage?.setItem(LAST_VIEW_STORAGE_KEY, view);
    } catch {
      // Ignore storage failures (private mode, disabled storage) — non-essential.
    }
  }, [view]);

  return (
    <main className="app-shell">
      <header className="app-header">
        <div className="app-identity">
          <div className="app-logo"><span>ES</span></div>
          <div><Text as="h1">Exercise Signature Explorer</Text><Text as="p">Track 2 · Omic Discordance Explained</Text></div>
        </div>
        <div className="science-disclaimer"><span>RESEARCH USE</span> Cross-cohort association evidence — not a treatment claim</div>
      </header>
      <nav className="view-nav" aria-label="Views">
        {NAV_VIEWS.map((item) => (
          <Button
            className={`view-nav__item ${view === item.id ? "view-nav__item--active" : ""}`}
            onClick={() => setView(item.id)}
            key={item.id}
          >
            <span>{item.eyebrow}</span>{item.label}
          </Button>
        ))}
      </nav>
      {view === "discordance" && <Discordance />}
      {view === "live" && <LiveDashboard />}
      {view === "generalized" && <GeneralizedQuery />}
      {view === "metabolomics" && <MetabolomicsCaseStudy />}
      {view === "methods" && <AboutMethods />}
      <footer className="app-footer"><span>System handoff · v0.1</span><span>Cross-cohort comparison · transparent assumptions · reproducible export</span></footer>
    </main>
  );
}
