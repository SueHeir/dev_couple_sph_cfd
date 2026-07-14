"""Independent failure-mode tests for the external-PSI eligibility firewall.

These tests deliberately do not use simulation output.  They demonstrate that a
plausible-looking manifest is still rejected when it omits a source table or a
measurement uncertainty, and that the current mismatched executable remains
ineligible.
"""
from __future__ import annotations

import json
import pathlib
import subprocess
import tempfile
import unittest


HERE = pathlib.Path(__file__).resolve().parent
CONTRACT = HERE / "validation_contract.py"
AUDIT = HERE.parent / "external_reference_audit.py"


class ValidationFirewallTests(unittest.TestCase):
    def test_current_case_fails_closed(self) -> None:
        result = subprocess.run(
            ["python3", str(AUDIT)], text=True, capture_output=True, check=False
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("EXTERNAL PSI COMPARISON: INELIGIBLE", result.stdout)

    def test_manifest_without_source_table_is_rejected(self) -> None:
        manifest = {
            "citation": "not evidence", "geometry": "circular_impinging_jet",
            "nozzle": "0.01 m", "stand_off_m": 0.08, "gas": "nitrogen",
            "material": "quartz", "grain_range_m": [0.0002, 0.0006],
            "forcing": "34 m/s", "duration_s": 100, "observable": "crater_depth",
            "units": "m", "uncertainty": {"absolute": 0.001},
            "data_file": "missing.csv", "adversarial_control": "wrong drag sign",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "external-comparison.json"
            path.write_text(json.dumps(manifest))
            result = subprocess.run(
                ["python3", str(CONTRACT), str(path)],
                text=True, capture_output=True, check=False,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("source data file is absent", result.stdout)


if __name__ == "__main__":
    unittest.main()
