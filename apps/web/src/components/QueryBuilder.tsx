import { createElement, useMemo, useState, type ReactNode } from "react";
import AvailabilityState from "./AvailabilityState";
import {
  getStudyContext,
  layerCapabilities,
  reconcileQueryForSpecies,
  resolveContext,
  type AnalysisQuery,
  type MotrpacSpecies,
  type OmicLayer,
  type SignatureSourceSpecies,
} from "../domain/analysis";

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

export const initialAnalysisQuery: AnalysisQuery = {
  signatureId: "pah-skm-v1.2",
  signatureSourceSpecies: "human",
  targetMotrpacSpecies: "rat",
  selectedOmics: ["transcriptomics", "proteomics"],
  tissue: "SKM-VL",
  sex: "female",
  timepoint: "8w",
  significanceThreshold: 0.05,
  includeNonSignificantResults: true,
};

function title(value: string) {
  return value.replaceAll("_", " ").replace(/\b\w/g, (character) => character.toUpperCase());
}

function Field({
  label,
  children,
}: {
  label: string;
  children: ReactNode;
}) {
  return (
    <div className="query-field">
      <span className="query-field__label">{label}</span>
      <div className="query-options">{children}</div>
    </div>
  );
}

export default function QueryBuilder({
  onQueryChange,
}: {
  onQueryChange: (query: AnalysisQuery) => void;
}) {
  const [query, setQuery] = useState(initialAnalysisQuery);
  const [changeNotice, setChangeNotice] = useState<string | null>(null);
  const context = getStudyContext(query.targetMotrpacSpecies)!;
  const availability = useMemo(
    () =>
      resolveContext({
        targetMotrpacSpecies: query.targetMotrpacSpecies,
        selectedOmics: query.selectedOmics,
        tissue: query.tissue,
        sex: query.sex,
        timepoint: query.timepoint,
      }),
    [query],
  );

  function update(next: AnalysisQuery) {
    setQuery(next);
    onQueryChange(next);
  }

  function changeSpecies(species: MotrpacSpecies) {
    if (species === query.targetMotrpacSpecies) return;
    const nextContext = getStudyContext(species)!;
    const previousContext = context;
    const reconciled = reconcileQueryForSpecies(query, species);
    setChangeNotice(
      `Context changed: ${previousContext.displayStudyDesign} → ${nextContext.displayStudyDesign}. Updated incompatible fields: ${reconciled.changedFields.join(", ")}.`,
    );
    update(reconciled.query);
  }

  function toggleOmic(layer: OmicLayer) {
    const exists = query.selectedOmics.includes(layer);
    const selectedOmics = exists
      ? query.selectedOmics.filter((item) => item !== layer)
      : [...query.selectedOmics, layer];
    update({ ...query, selectedOmics });
  }

  const structurallyInvalid = query.selectedOmics.length === 0;

  return (
    <div className="query-builder">
      <div className="query-builder__title">
        <div>
          <span className="query-builder__eyebrow">Analysis query</span>
          <strong>User selections</strong>
        </div>
        <span className="query-builder__version">CATALOG 1.0</span>
      </div>

      <Field label="Disease-signature source species">
        {(["human", "rat", "mouse", "other"] as SignatureSourceSpecies[]).map((species) => (
          <Button
            className={`query-chip ${query.signatureSourceSpecies === species ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, signatureSourceSpecies: species })}
            key={species}
          >
            {title(species)}
          </Button>
        ))}
      </Field>

      <Field label="Target MoTrPAC species">
        {(["human", "rat"] as MotrpacSpecies[]).map((species) => {
          const speciesContext = getStudyContext(species)!;
          return (
            <Button
              className={`species-choice ${query.targetMotrpacSpecies === species ? "species-choice--selected" : ""}`}
              onClick={() => changeSpecies(species)}
              key={species}
            >
              <strong>{title(species)}</strong>
              <span>{speciesContext.displayStudyDesign}</span>
            </Button>
          );
        })}
      </Field>

      <div className="derived-study-context">
        <span>Derived study design</span>
        <strong>{context.displayStudyDesign}</strong>
        <small>Read-only catalog metadata · not a user input</small>
      </div>

      {changeNotice && <div className="query-change-notice">{changeNotice}</div>}

      <Field label="Omic layers · multi-select">
        {layerCapabilities.map((capability) => (
          <Button
            className={`omic-choice ${query.selectedOmics.includes(capability.layer) ? "omic-choice--selected" : ""} ${capability.status !== "implemented" ? "omic-choice--planned" : ""}`}
            onClick={() => toggleOmic(capability.layer)}
            key={capability.layer}
          >
            <span>{title(capability.layer)}</span>
            <b>{capability.status === "implemented" ? "Implemented" : "Planned"}</b>
          </Button>
        ))}
      </Field>

      <Field label="Tissue">
        {context.tissues.map((tissue) => (
          <Button className={`query-chip ${query.tissue === tissue ? "query-chip--selected" : ""}`} onClick={() => update({ ...query, tissue })} key={tissue}>
            {tissue}
          </Button>
        ))}
      </Field>

      <Field label="Sex">
        {context.sexes.map((sex) => (
          <Button className={`query-chip ${query.sex === sex ? "query-chip--selected" : ""}`} onClick={() => update({ ...query, sex })} key={sex}>
            {title(sex)}
          </Button>
        ))}
      </Field>

      <Field label="Time point">
        {context.timepoints.map((timepoint) => (
          <Button className={`query-chip ${query.timepoint === timepoint ? "query-chip--selected" : ""}`} onClick={() => update({ ...query, timepoint })} key={timepoint}>
            {title(timepoint)}
          </Button>
        ))}
      </Field>

      <Field label="Adjusted-significance threshold">
        {[0.01, 0.05, 0.1].map((threshold) => (
          <Button
            className={`query-chip ${query.significanceThreshold === threshold ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, significanceThreshold: threshold })}
            key={threshold}
          >
            q &lt; {threshold}
          </Button>
        ))}
      </Field>

      <div className="query-toggle-row">
        <div><strong>Include nonsignificant measurements</strong><span>Preserve effects and adjusted p-values</span></div>
        <Button
          className={`query-toggle ${query.includeNonSignificantResults ? "query-toggle--on" : ""}`}
          onClick={() => update({ ...query, includeNonSignificantResults: !query.includeNonSignificantResults })}
          aria-pressed={query.includeNonSignificantResults}
        >
          <span />
        </Button>
      </div>

      <AvailabilityState result={availability} />

      <Button className="run-comparison" disabled={structurallyInvalid}>
        {structurallyInvalid ? "Select at least one omic layer" : "Run comparison →"}
      </Button>
      <div className="run-note">Availability never disables submission. Only a structurally invalid query does.</div>
    </div>
  );
}
