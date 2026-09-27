"""Robustness matrix, metabolite rows: RefMet names, HMDB ids, unknown names, through the CLI."""
import subprocess
import sys

import pandas as pd
import pytest

from motrpac_probe import render

CASES = {
    "refmet_names": "refmet_name,direction\nCitric acid,-1\nSuccinic acid,-1\nFumaric acid,-1\nMalic acid,-1\nLactic acid,1\n",
    "hmdb_ids": "hmdb,direction\nHMDB0000094,-1\nHMDB0000254,-1\nHMDB0000134,-1\nHMDB0000156,-1\nHMDB0000190,1\n",
    "unknown_names": "refmet_name,direction\nCitric acid,-1\nSuccinic acid,-1\nNotAMetabolite,1\nFumaric acid,-1\n",
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_metab_signature_cli(tmp_path, case):
    f = tmp_path / f"{case}.csv"
    f.write_text(CASES[case])
    out = tmp_path / "out"
    p = subprocess.run([sys.executable, "-m", "motrpac_probe.cli", "run", "--quiet", "--nboot", "200",
                        "--signature", str(f), "--out", str(out)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]
    assert render.validate_report(out / "report.html") == []
    mp = pd.read_csv(out / "tables" / "mapping.csv", dtype=str)
    sc = pd.read_csv(out / "tables" / "column_scores.csv")
    assert sc.column_id.str.contains("\\|METAB\\|").all()
    if case == "unknown_names":
        assert (mp.status != "mapped").sum() == 1
    else:
        assert (mp.status == "mapped").sum() >= 4
