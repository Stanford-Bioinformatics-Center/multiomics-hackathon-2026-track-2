"""Request/response contract for the analysis service.

Plain dataclasses so the service has no web-framework dependency. Gate 5 generates Pydantic
models + an OpenAPI spec from these (and converges toward the React app's domain/analysis.ts
names). Everything here is JSON-safe: only str / int / float / bool / None / list / dict.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SignatureRow:
    """One input row of a directed disease signature (as uploaded / pasted / from a fixture)."""
    gene_symbol: Optional[str] = None
    uniprot: Optional[str] = None
    ensembl: Optional[str] = None
    rat_symbol: Optional[str] = None
    direction: Optional[str] = None   # "+1"/"-1"/"up"/"down"; parsed by the engine
    group: Optional[str] = None
    weight: Optional[str] = None       # carried through; NOT used as a statistical weight
    source: Optional[str] = None


@dataclass
class AnalysisRequest:
    """A researcher's analysis request. The exercise study design is DERIVED from
    target_species (rat -> chronic training, human -> acute bout); there is no acute/chronic input.
    'include_nonsignificant' is presentation-only and MUST NOT change the run (not in run_id)."""
    # signature: exactly one of signature_rows / signature_csv_text / example_name
    signature_rows: Optional[list[SignatureRow]] = None
    signature_csv_text: Optional[str] = None
    example_name: Optional[str] = None          # a built-in example under hackathon/tool/examples
    signature_name: Optional[str] = None

    target_species: str = "rat"                  # "rat" | "human"
    selected_omics: list[str] = field(default_factory=lambda: ["transcriptomics", "proteomics"])
    tissue: Optional[str] = None                 # display filter only
    sex: Optional[str] = None                    # display filter only
    timepoint: Optional[str] = None              # display filter only
    fdr_threshold: float = 0.05                  # analysis parameter -> part of run_id
    include_nonsignificant: bool = True          # presentation-only -> NOT part of run_id


# ---- response ------------------------------------------------------------------------------------

@dataclass
class MappingAuditRow:
    input_row: int
    input_id: str
    selected_symbol: str
    selected_method: str
    status: str
    n_candidates: int
    ambiguous: bool
    candidates: str


@dataclass
class ContrastPair:
    """Both contrasts are first-class, never reduced to 'healthy gene set'."""
    disease_numerator: str
    disease_denominator: str
    exercise_numerator: str
    exercise_denominator: str


@dataclass
class FeatureEvidence:
    """One signature gene in one comparison column. Every value traces to a store record so that
    point -> evidence row -> source record is reconstructable (lineage fields below)."""
    evidence_id: str
    gene: str
    disease_direction: int
    column_id: str
    comparison_label: str
    dataset: str
    species: str
    tissue: str
    layer: str
    sex: str
    timepoint: str
    exercise_logfc: Optional[float]
    exercise_stat: Optional[float]
    fdr_bh: Optional[float]          # within-comparison feature-level BH (from the store)
    opposed: Optional[bool]          # sign(exercise logFC) opposes disease direction
    measured: bool
    # ---- lineage (point -> evidence -> source) ----
    mapping_decision_id: str = ""    # ties this gene back to a mapping_audit row
    source_feature_id: str = ""      # the store feature_id backing this gene in this column
    n_collapsed: Optional[int] = None  # how many assay features collapsed into this gene
    aggregation_method: str = ""     # e.g. "max_abs_stat" when n_collapsed > 1, else "single_feature"


@dataclass
class ColumnResult:
    """Set-level result for one comparison column, under the frozen multiplicity family."""
    column_id: str
    comparison_label: str
    dataset: str
    species: str
    tissue: str
    layer: str
    sex: str
    timepoint: str
    n_measured: int
    n_opposed: int
    n_same: int
    camera_t: Optional[float]
    camera_p: Optional[float]         # unadjusted
    camera_fdr: Optional[float]       # BH across the whole family (NOT re-run on filtering)
    camera_bonferroni: Optional[float] # Bonferroni across that same pre-filter family
    verdict: str                      # "opposed" | "same direction" | "no set-level shift"
    significant: bool                 # camera_fdr < fdr_threshold


@dataclass
class LayerDiscordanceRow:
    comparison_label: str
    tissue: str
    genes_in_both: int
    genome_wide_rho: Optional[float]
    signature_rho: Optional[float]
    rna_opposed: int
    protein_opposed: int
    both_opposed: int
    concordant: int
    rna_only: int
    protein_only: int
    opposite: int
    neither: int


@dataclass
class MultiplicityFamily:
    method: str                       # "BH"
    methods: list[str]                # BH and Bonferroni, on the same primary cameraPR family
    threshold: float
    family_size: int                  # number of columns in the family (e.g. 52)
    n_tests: int                      # non-missing cameraPR tests actually adjusted
    column_ids: list[str]             # all ordered column IDs in the family


@dataclass
class AnalysisResponse:
    run_id: str
    schema_version: str
    status: str                       # "ok" | "no_compatible_data" | "no_mapped_features"
    message: str
    signature_name: str
    n_input_rows: int
    n_counted_genes: int
    n_mapped_rows: int
    returned_layers: list[str]        # layers actually present in the results (drives viz mode)
    contrasts: dict                   # ContrastPair per dataset, as plain dict
    mapping_audit: list[MappingAuditRow]
    columns: list[ColumnResult]
    features: list[FeatureEvidence]
    layer_discordance: list[LayerDiscordanceRow]
    multiplicity_family: MultiplicityFamily
    headline: dict                    # the frozen headline result (male rat SKM-GN protein 8 wk, etc.)
    caveats: list[str]
    guardrails: list[dict]            # Safe/Unsafe rows
    provenance: dict
