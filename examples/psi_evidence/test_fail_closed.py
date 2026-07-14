#!/usr/bin/env python3
"""Regression test for the claim guard, not scientific validation."""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent


class PsiEvidenceGuardTests(unittest.TestCase):
    def test_absent_external_package_blocks_claim(self) -> None:
        result = subprocess.run(
            [sys.executable, str(HERE / "audit.py")],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("PSI EVIDENCE: BLOCKED", result.stdout)


if __name__ == "__main__":
    unittest.main()
