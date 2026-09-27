"""Family-wise correction is computed by the engine before display filtering."""

import numpy as np
import pytest

from motrpac_probe import core


def test_camera_bonferroni_uses_same_nonmissing_family_as_bh(monkeypatch):
    p_by_id = {"a": 0.001, "b": 0.03, "c": np.nan, "d": 0.8}

    def score_column(_store, cid, _sig, **_kw):
        return {"column_id": cid, "camera_p": p_by_id[cid]}

    monkeypatch.setattr(core, "score_column", score_column)
    scored = core.score_columns(None, list(p_by_id), None).set_index("column_id")

    assert scored.loc["a", "camera_bonferroni"] == pytest.approx(0.003)
    assert scored.loc["b", "camera_bonferroni"] == pytest.approx(0.09)
    assert np.isnan(scored.loc["c", "camera_bonferroni"])
    assert scored.loc["d", "camera_bonferroni"] == pytest.approx(1.0)
    assert scored.loc["b", "camera_fdr"] == pytest.approx(0.045)

    # A display subset must retain the precomputed correction, not use two rows as the family.
    displayed = scored.loc[["a", "b"]]
    assert displayed.loc["b", "camera_bonferroni"] == pytest.approx(0.09)
