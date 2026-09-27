import type { ReactNode } from "react";

type NodeProps = {
  x: number;
  y: number;
  width: number;
  height?: number;
  title: string;
  detail?: string;
  tone?: "navy" | "purple" | "orange" | "teal" | "gray" | "coral";
  dashed?: boolean;
};

function Node({
  x,
  y,
  width,
  height = 74,
  title,
  detail,
  tone = "navy",
  dashed = false,
}: NodeProps) {
  return (
    <g className={`flow-node flow-node--${tone} ${dashed ? "flow-node--dashed" : ""}`}>
      <rect x={x} y={y} width={width} height={height} rx="16" />
      <foreignObject x={x + 12} y={y + 8} width={width - 24} height={height - 16}>
        <div className="flow-node__content">
          <strong>{title}</strong>
          {detail && <span>{detail}</span>}
        </div>
      </foreignObject>
    </g>
  );
}

function Edge({
  d,
  label,
  x,
  y,
  width = 200,
  planned = false,
}: {
  d: string;
  label: string;
  x: number;
  y: number;
  width?: number;
  planned?: boolean;
}) {
  return (
    <g className={`flow-connector ${planned ? "flow-connector--planned" : ""}`}>
      <path d={d} markerEnd={planned ? "url(#technical-arrow-planned)" : "url(#technical-arrow)"} />
      <foreignObject x={x - width / 2} y={y - 13} width={width} height="30">
        <div className="flow-connector__label">{label}</div>
      </foreignObject>
    </g>
  );
}

function LegendItem({ className, children }: { className: string; children: ReactNode }) {
  return <span><i className={className} />{children}</span>;
}

export default function TechnicalFlowDiagram() {
  return (
    <>
      <div className="flow-legend">
        <LegendItem className="flow-key flow-key--solid">Implemented in MVP</LegendItem>
        <LegendItem className="flow-key flow-key--planned">Planned adapter or relationship</LegendItem>
        <LegendItem className="flow-key flow-key--human">Researcher confirmation or scientific empty state</LegendItem>
        <span>Every arrow names its payload or interaction</span>
      </div>

      <div className="technical-flow-shell">
        <svg className="technical-flow technical-flow--expanded" viewBox="0 0 1600 1180" role="img" aria-labelledby="technical-flow-title technical-flow-description">
          <title id="technical-flow-title">Exercise Signature Explorer technical architecture and data flow</title>
          <desc id="technical-flow-description">
            A researcher submits a directed disease signature, its source species, a target MoTrPAC species, and one or more omic layers. The catalog derives study design from species. Query validation and dataset availability precede identifier mapping. Compatible contexts continue through independent omic adapters into multi-omic evidence integration. Incompatible contexts produce a scientific empty state and alternatives. The optional language model remains outside the calculation path.
          </desc>
          <defs>
            <marker id="technical-arrow" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
              <path d="M0,0 L0,6 L9,3 z" className="flow-arrowhead" />
            </marker>
            <marker id="technical-arrow-planned" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto" markerUnits="strokeWidth">
              <path d="M0,0 L0,6 L9,3 z" className="flow-arrowhead flow-arrowhead--planned" />
            </marker>
          </defs>

          <text x="36" y="34" className="flow-column-label">USER + QUERY CONTRACT</text>
          <text x="410" y="34" className="flow-column-label">CATALOG + AVAILABILITY</text>
          <text x="805" y="34" className="flow-column-label">SCIENTIFIC PROCESSING</text>
          <text x="1285" y="34" className="flow-column-label">EVIDENCE DELIVERY</text>

          <Node x={36} y={65} width={290} title="Researcher" detail="confirms every structured selection" />
          <Node x={36} y={185} width={290} title="Directed-signature input" detail="directed features + source metadata + source species" tone="orange" />
          <Node x={36} y={305} width={290} title="MoTrPAC species selector" detail="Human or Rat · study design is not editable" />
          <Node x={36} y={425} width={290} title="Omic multi-select" detail="one or more layers · implemented status visible" />
          <Node x={36} y={545} width={290} title="Context filters" detail="tissue · sex · time point · q threshold" />
          <Node x={36} y={665} width={290} title="Query validation" detail="schema + required fields + explicit confirmation" />

          <Node x={430} y={65} width={310} title="Study-context resolver" detail="Human → acute · Rat → chronic" tone="purple" />
          <Node x={430} y={185} width={310} title="Dataset / context catalog" detail="species · design · tissues · sexes · time points · versions" tone="orange" />
          <Node x={430} y={305} width={310} title="Layer capability registry" detail="implemented · partial · planned · entity types" tone="orange" />
          <Node x={430} y={545} width={310} title="Dataset Availability Resolver" detail="deterministic catalog lookup before analysis" tone="purple" />
          <Node x={430} y={740} width={310} title="Scientific empty state" detail="unsupported layer or no matching tissue / time point" tone="coral" />
          <Node x={430} y={860} width={310} title="Available alternatives" detail="supported layers and compatible catalog filters" tone="coral" />

          <Node x={825} y={65} width={310} title="Identifier / ortholog normalization" detail="aliases · namespaces · canonical IDs · human ↔ rat" tone="purple" />
          <Node x={825} y={185} width={310} title="Mapping confirmation gate" detail="recognized · remapped · ambiguous · unmapped" tone="coral" />

          <Node x={795} y={345} width={175} title="Genomics adapter" detail="variant / gene" tone="gray" dashed />
          <Node x={990} y={345} width={175} title="Epigenomics adapter" detail="regulatory context" tone="gray" dashed />
          <Node x={795} y={455} width={175} title="Transcriptomics adapter" detail="gene / transcript" tone="purple" />
          <Node x={990} y={455} width={175} title="Proteomics adapter" detail="protein / PTM site" tone="purple" />
          <Node x={893} y={565} width={175} title="Metabolomics adapter" detail="metabolite / pathway" tone="gray" dashed />

          <Node x={825} y={705} width={340} title="Multi-Omic Evidence Integration" detail="independent layers · pairwise relations · missingness retained" tone="purple" />
          <Node x={825} y={825} width={340} title="Comparison + classification" detail="disease opposition + eight RNA/protein states" tone="purple" />
          <Node x={825} y={945} width={340} title="Versioned evidence assembly" detail="effects · q-values · classifications · provenance · limits" tone="teal" />

          <Node x={1280} y={65} width={280} title="Optional LLM" detail="outside scientific calculation path" tone="gray" dashed />
          <Node x={1280} y={185} width={280} title="Query proposal / explanation" detail="researcher must confirm before execution" tone="gray" dashed />
          <Node x={1280} y={545} width={280} title="Results + provenance API" detail="features · methods · mappings · limitations" tone="teal" />
          <Node x={1280} y={665} width={280} title="Interactive dashboard" detail="single-layer · pairwise · heatmap / matrix" tone="teal" />
          <Node x={1280} y={785} width={280} title="Scientific result states" detail="mapped-none · weak evidence · measured missingness" tone="teal" />
          <Node x={1280} y={945} width={280} title="Evidence-bundle export" detail="CSV + JSON + manifest + provenance" tone="teal" />

          <Edge d="M181 139 V185" label="directed features + source metadata" x={181} y={162} width={230} />
          <Edge d="M326 342 H430 V102" label="target MoTrPAC species" x={378} y={286} />
          <Edge d="M585 139 V185" label="resolved species request" x={585} y={162} />
          <Edge d="M585 259 V545" label="study availability and metadata" x={585} y={402} width={225} />
          <Edge d="M430 102 H380 V462 H326" label="derived acute / chronic context" x={380} y={221} width={215} />
          <Edge d="M430 342 H385 V462 H326" label="layer statuses" x={385} y={364} width={130} />
          <Edge d="M181 499 V545" label="selected omics" x={181} y={522} />
          <Edge d="M181 619 V665" label="tissue + sex + time point" x={181} y={642} width={210} />
          <Edge d="M326 702 H375 V582 H430" label="availability lookup" x={378} y={650} />
          <Edge d="M585 619 V740" label="no-match reason + alternatives" x={585} y={679} width={225} />
          <Edge d="M585 814 V860" label="compatible catalog options" x={585} y={837} />
          <Edge d="M430 897 H365 V342 H326" label="researcher revises query" x={365} y={820} width={185} />

          <Edge d="M740 582 H780 V102 H825" label="validated query context" x={780} y={288} width={200} />
          <Edge d="M980 139 V185" label="canonical + ortholog identifiers" x={980} y={162} width={225} />
          <Edge d="M980 259 V300 H882 V345" label="confirmed IDs by entity type" x={882} y={301} width={200} />
          <Edge d="M980 259 V300 H1077 V345" label="planned relationships" x={1077} y={301} width={175} planned />

          <Edge d="M882 419 V705" label="planned variant relationships" x={830} y={554} width={190} planned />
          <Edge d="M1077 419 V705" label="planned regulatory relationships" x={1127} y={554} width={205} planned />
          <Edge d="M882 529 V705" label="matched transcript effects" x={852} y={616} width={190} />
          <Edge d="M1077 529 V705" label="matched protein effects" x={1112} y={616} width={180} />
          <Edge d="M980 639 V705" label="planned metabolite relationships" x={980} y={672} width={210} planned />
          <Edge d="M995 779 V825" label="matched omic effects" x={995} y={802} />
          <Edge d="M995 899 V945" label="classification results" x={995} y={922} />

          <Edge d="M1165 982 H1220 V582 H1280" label="versioned evidence bundle" x={1220} y={906} width={195} />
          <Edge d="M1420 619 V665" label="deterministic results" x={1420} y={642} />
          <Edge d="M1420 739 V785" label="effects + q-values + missingness" x={1420} y={762} width={220} />
          <Edge d="M1165 982 H1280" label="versioned evidence bundle" x={1222} y={960} width={190} />

          <Edge d="M1420 139 V185" label="user intent only" x={1420} y={162} planned />
          <Edge d="M1280 222 H1215 V102 H1135" label="proposed structured query" x={1215} y={154} width={185} planned />
          <Edge d="M1280 222 H1200 V702 H1280" label="explain empty state or validated result" x={1200} y={450} width={225} planned />
        </svg>
      </div>

      <div className="technical-description">
        <strong>Accessible flow summary.</strong>
        <p>
          The researcher chooses the disease-signature source species and the target MoTrPAC species separately. The target species resolves a read-only study design from the versioned catalog: the current human context is acute and the current rat context is chronic training. Availability validation occurs before mapping. Compatible queries continue through canonical identifier and ortholog normalization, then independent transcriptomic and proteomic adapters. Planned genomics, epigenomics, and metabolomics adapters use dashed styling and do not represent functioning analysis paths. A no-match returns a scientific empty state with alternatives, not a technical error.
        </p>
        <p>
          Independent entity types include variants, genes, transcripts, proteins, PTM sites, and metabolites. The MVP validates gene/transcript-to-protein relationships only. Other regulatory and pathway relationships remain planned, and no direct one-to-one match is implied across every layer.
        </p>
      </div>
    </>
  );
}
