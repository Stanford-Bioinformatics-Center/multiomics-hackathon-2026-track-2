"""Offline tests for the accession-only Workbench fetcher."""

from __future__ import annotations

import hashlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fetch_workbench import download_study, validate_accession  # noqa: E402


class FakeResponse(io.BytesIO):
    def __init__(self, payload: bytes, url: str):
        super().__init__(payload)
        self.url = url

    def geturl(self) -> str:
        return self.url


class FetchWorkbenchTests(unittest.TestCase):
    def test_accession_only(self) -> None:
        self.assertEqual(validate_accession("st000763"), "ST000763")
        for invalid in ("ST763", "https://example.org", "ST123456/../data", "ST123456?x=1"):
            with self.assertRaises(ValueError):
                validate_accession(invalid)

    def test_download_writes_official_payloads_and_hashes(self) -> None:
        factors = json.dumps([{"study_id": "ST123456", "local_sample_id": "S1",
                               "factors": "Group: Case"}]).encode()
        data = json.dumps([{"study_id": "ST123456", "analysis_id": "AN1",
                            "metabolite_id": "ME1", "DATA": {"S1": 1}}]).encode()
        payloads = {"factors": factors, "data": data}

        def fake_urlopen(url: str, timeout: int) -> FakeResponse:
            self.assertEqual(timeout, 60)
            endpoint = url.rsplit("/", 1)[-1]
            return FakeResponse(payloads[endpoint], url)

        with tempfile.TemporaryDirectory() as temp, patch("fetch_workbench.urlopen", side_effect=fake_urlopen):
            target = download_study("ST123456", Path(temp))
            manifest = json.loads((target / "source_manifest.json").read_text())
            for name, payload in payloads.items():
                self.assertEqual((target / f"{name}.json").read_bytes(), payload)
                self.assertEqual(manifest["files"][name]["sha256"], hashlib.sha256(payload).hexdigest())
            with self.assertRaises(FileExistsError):
                download_study("ST123456", Path(temp))

    def test_redirect_to_other_host_is_rejected(self) -> None:
        payload = json.dumps([{"study_id": "ST123456"}]).encode()
        with tempfile.TemporaryDirectory() as temp, patch(
            "fetch_workbench.urlopen",
            return_value=FakeResponse(payload, "https://example.org/redirect"),
        ):
            with self.assertRaisesRegex(ValueError, "Unexpected redirect"):
                download_study("ST123456", Path(temp))
            self.assertFalse((Path(temp) / "ST123456" / "factors.json").exists())


if __name__ == "__main__":
    unittest.main()
