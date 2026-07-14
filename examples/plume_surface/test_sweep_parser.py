#!/usr/bin/env python3
"""Contract test for parsing the executable's current regression report."""

import importlib.util
import pathlib
import unittest


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("sweep", HERE / "sweep.py")
assert SPEC and SPEC.loader
sweep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sweep)


CURRENT_REPORT = """
# SPH continuum (seam, MacDonald, drag+gradP) : 0.18310   <- the coupled-continuum measurement
# Wen & Yu 1966 correlation (packed-bed comparator) : 0.18850   rel.err 2.87%  (regression limit 20.0%)
# DEM discrete  (same seam, MacDonald, drag+gradP) : 0.18310   consistency rel.err 0.00%  (regression limit 10.0%)
# fault controls: omit-gradP 0.4100 (+117.8%)  eps-power-bug 0.0773 (-59.0%)  => both FAIL as required (must exceed regression limit 20.0%)
    0.50      0.0916     +1.0000e-04      3.0000e+01       0.3000   packed
    1.20      0.2197     +2.0000e-03      1.0000e+00       0.0100   FLUIDIZES
REGRESSION: PASS  (packed-bed force-transfer checks; not plume/crater validation)
"""


class SweepParserTests(unittest.TestCase):
    def test_current_report_vocabulary_is_parseable(self):
        data = sweep.parse_output(CURRENT_REPORT)
        self.assertEqual(data["result"].split()[1], "PASS")
        self.assertEqual(data["tol_wenyu"], 20.0)
        self.assertEqual(data["tol_dem"], 10.0)
        self.assertEqual(len(data["sweep"]), 2)


if __name__ == "__main__":
    unittest.main()
