#!/usr/bin/env python3
"""Fail closed for the pending PSI validation claim.

This script deliberately cannot certify a crater/erosion/ejecta validation.
File names, hashes, and self-declared manifest fields establish neither the
provenance of external data nor that an advancing CFD case matches it. Those
facts require a reviewed, independently runnable comparator. Consequently this
gate always exits nonzero; it guards against an unsupported claim, not a
numerical acceptance test.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


REQUIRED = {
    "citation", "source_url", "source_page_or_figure", "source_data_file",
    "source_data_sha256", "source_uncertainty", "geometry", "nozzle_d_m",
    "stand_off_m", "gas", "material", "grain_range_m", "duration_s",
    "observable", "units", "prediction_file", "prediction_sha256",
    "wrong_coupling_file", "wrong_coupling_sha256", "comparison_script",
}
ALLOWED_OBSERVABLES = {"crater_depth", "crater_volume", "ejecta_mass_rate"}


def blocked(message: str) -> None:
    print(f"PSI EVIDENCE: BLOCKED — {message}")
    raise SystemExit(1)


def digest(path: Path, expected: object, label: str) -> None:
    if not path.is_file():
        blocked(f"{label} is absent: {path.name}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        blocked(f"{label} checksum does not match its manifest entry")


def has_rows(path: Path, label: str) -> None:
    try:
        rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    except (csv.Error, UnicodeDecodeError) as exc:
        blocked(f"{label} is not readable CSV: {exc}")
    if not rows:
        blocked(f"{label} has no quantitative rows")


def main() -> None:
    here = Path(__file__).resolve().parent
    manifest_path = here / "external_comparison.json"
    if not manifest_path.is_file():
        blocked("no primary-source comparison manifest is committed")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        blocked(f"comparison manifest is invalid JSON: {exc}")
    missing = sorted(REQUIRED - manifest.keys())
    if missing:
        blocked("manifest lacks " + ", ".join(missing))
    if manifest["geometry"] != "circular_impinging_jet":
        blocked("reference is not a circular impinging-jet experiment")
    if manifest["observable"] not in ALLOWED_OBSERVABLES:
        blocked("observable is not a quantitative crater or ejecta metric")
    if not isinstance(manifest["source_uncertainty"], (int, float)) or manifest["source_uncertainty"] <= 0:
        blocked("source uncertainty must be a positive measured uncertainty")
    source = here / manifest["source_data_file"]
    prediction = here / manifest["prediction_file"]
    wrong = here / manifest["wrong_coupling_file"]
    comparator = here / manifest["comparison_script"]
    digest(source, manifest["source_data_sha256"], "primary-source data")
    digest(prediction, manifest["prediction_sha256"], "model prediction")
    digest(wrong, manifest["wrong_coupling_sha256"], "wrong-coupling control")
    has_rows(source, "primary-source data")
    has_rows(prediction, "model prediction")
    has_rows(wrong, "wrong-coupling control")
    if not comparator.is_file():
        blocked("independent comparator script is absent")
    blocked(
        "a manifest-shaped package is not independently verified external PSI evidence; "
        "do not claim validation until a reviewed comparator establishes matched inputs, "
        "uncertainty, and a failing wrong-coupling control"
    )


if __name__ == "__main__":
    main()
