import { createElement, type ReactNode } from "react";
import type { AvailabilityResult, OmicLayer } from "../domain/analysis";

function Text({
  as = "div",
  className = "",
  children,
}: {
  as?: "div" | "span" | "p" | "h3";
  className?: string;
  children: ReactNode;
}) {
  return createElement(as, { className }, children);
}

function formatLayer(layer: OmicLayer) {
  return layer.charAt(0).toUpperCase() + layer.slice(1);
}

export default function AvailabilityState({ result }: { result: AvailabilityResult }) {
  const context = result.resolvedContext;
  const statusLabel = {
    available: "Compatible context",
    partially_available: "Partial availability",
    unsupported_omics: "Requested layer unavailable",
    no_matching_context: "No matching context",
  }[result.status];

  return (
    <div className={`availability-state availability-state--${result.status}`}>
      <div className="availability-state__head">
        <span className="availability-state__mark">
          {result.status === "available" ? "✓" : result.status === "partially_available" ? "◐" : "—"}
        </span>
        <div>
          <Text className="availability-state__kicker">Deterministic catalog lookup</Text>
          <Text as="h3">{statusLabel}</Text>
        </div>
      </div>
      {context && (
        <div className="resolved-context">
          <span>{context.species}</span>
          <span>{context.displayStudyDesign}</span>
          <span>{context.datasetVersion}</span>
        </div>
      )}
      <Text as="p">{result.message}</Text>
      <div className="availability-layers">
        {result.availableOmics.map((layer) => (
          <span className="availability-layer availability-layer--yes" key={layer}>
            ✓ {formatLayer(layer)}
          </span>
        ))}
        {result.unavailableOmics.map((layer) => (
          <span className="availability-layer availability-layer--no" key={layer}>
            — {formatLayer(layer)} · planned
          </span>
        ))}
      </div>
      {result.status === "no_matching_context" && context && (
        <div className="availability-actions">
          Available: {context.tissues.join(" / ")} · {context.timepoints.join(" / ")}
        </div>
      )}
    </div>
  );
}
