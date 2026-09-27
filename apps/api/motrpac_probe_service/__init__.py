"""UI-independent service adapter over the motrpac_probe scientific engine.

The single entry point is `run_analysis(request) -> AnalysisResponse`. It performs no statistics;
every number comes from `motrpac_probe`. See DECISIONS.md ADR-0001/0004/0007 and apps/api/README.md.
"""
from .schema import AnalysisRequest, AnalysisResponse, SignatureRow
from .service import SCHEMA_VERSION, run_analysis

__all__ = ["AnalysisRequest", "AnalysisResponse", "SignatureRow", "run_analysis", "SCHEMA_VERSION"]
