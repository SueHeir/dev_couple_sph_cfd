#!/usr/bin/env python3
"""Reject PSI validation claims that lack a comparable external data manifest.

This is deliberately independent of simulation output.  A comparison is only
eligible after its source data declares the experimental regime and uncertainty
needed to decide whether the two quantities are commensurate.  It prevents a
plotting driver from converting a convenient paper citation into a validation
gate.
"""
from __future__ import annotations

import json
import pathlib
import sys


REQUIRED = {
    "citation", "geometry", "nozzle", "stand_off_m", "gas", "material",
    "grain_range_m", "forcing", "duration_s", "observable", "units",
    "uncertainty", "data_file", "adversarial_control",
}
SAME_OBSERVABLE = {"crater_depth", "crater_volume", "ejecta_mass_rate"}


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: validation_contract.py <external-comparison.json>")
    manifest_path = pathlib.Path(sys.argv[1])
    if not manifest_path.is_file():
        print("EXTERNAL PSI COMPARISON: INELIGIBLE")
        print("INELIGIBLE: no external-comparison manifest with source data is committed")
        raise SystemExit(1)
    manifest = json.loads(manifest_path.read_text())
    missing = sorted(REQUIRED - manifest.keys())
    if missing:
        print("EXTERNAL PSI COMPARISON: INELIGIBLE")
        print("INELIGIBLE: manifest is missing " + ", ".join(missing))
        raise SystemExit(1)
    data_path = manifest_path.parent / manifest["data_file"]
    failures = []
    if manifest["geometry"] != "circular_impinging_jet":
        failures.append("reference geometry is not a circular impinging jet")
    if manifest["observable"] not in SAME_OBSERVABLE:
        failures.append("observable is not crater depth, volume, or ejecta mass rate")
    if not isinstance(manifest["uncertainty"], dict) or "absolute" not in manifest["uncertainty"]:
        failures.append("reference uncertainty is not a declared absolute measurement uncertainty")
    if not data_path.is_file():
        failures.append(f"source data file is absent: {data_path.name}")
    if failures:
        print("EXTERNAL PSI COMPARISON: INELIGIBLE")
        for failure in failures:
            print("INELIGIBLE:", failure)
        raise SystemExit(1)
    print("EXTERNAL PSI COMPARISON: ELIGIBLE")


if __name__ == "__main__":
    main()
