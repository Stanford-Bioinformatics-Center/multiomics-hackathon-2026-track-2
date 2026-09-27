"""Fixed caveat text and the Safe / Unsafe guardrails table (printed in every report and in README)."""

# From the team's planning workbook (sheet "Start Here", section 4 "Interpretation guardrails"; file named in
# README), with the disease generalised to "the disease signature", plus two rows from our own results.
GUARDRAILS = [
    ("A disease-signature pathway is directionally opposed after exercise.", "Exercise reverses or treats the disease."),
    ("Blood and muscle share a pathway label.", "The same genes changed or one tissue signalled to the other."),
    ("An early phosphosite and later target-motif signal are temporally compatible.",
     "The phosphosite caused the later transcriptional response."),
    ("A plasma protein matches tissue RNA/protein in time and direction.",
     "The tissue is proven to be the source of the plasma protein."),
    ("A PLIER LV is associated with a pathway prior.", "The pathway caused the LV or independently validates it."),
    ("RNA for a gene rose after one bout.", "The protein rose, or will rise."),
    ('A gene is "protein-only" in one contrast.', "The gene is post-transcriptionally regulated."),
]


def fixed_caveats():
    return [
        "'Opposed' is a sign comparison between a disease-versus-healthy difference and an exercise-versus-control "
        "difference in different people or animals. It is not a therapeutic claim.",
        "All MoTrPAC statistics here are group-level contrasts across genes, not per-person or per-animal "
        "measurements: we cannot tell whether the individuals whose RNA rose are those whose protein rose.",
        "Cross-species and cross-tissue: rat genes are mapped to human orthologs (RGD); rat protein is gastrocnemius "
        "while human biopsies are vastus lateralis; rat samples were taken about 48 h after the last training bout "
        "(adaptation, not a bout response).",
        "A random gene set from the same pathway class can look just as opposed (section D). Read the class-null "
        "percentile before attributing the result to the disease signature itself.",
    ]
