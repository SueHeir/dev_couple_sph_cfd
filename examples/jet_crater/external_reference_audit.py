#!/usr/bin/env python3
"""Fail closed when an external PSI trace is not comparable to this case.

This is an evidence audit, not a validator.  Exit status 1 is the correct
result while the executable and the published experiment have different
geometry/material/time scales.  Keeping it executable prevents an exploratory
plot from silently being presented as a crater-data comparison.
"""
from __future__ import annotations

import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
CFG = (HERE / "config.toml").read_text()


def value(key: str) -> float:
    match = re.search(rf"^{re.escape(key)}\s*=\s*([0-9.eE+-]+)", CFG, re.MULTILINE)
    if not match:
        raise RuntimeError(f"missing {key} in config.toml")
    return float(match.group(1))


def main() -> None:
    # The generic manifest gate is intentionally evaluated first.  This branch
    # has no source-data manifest because the only retrieved paper is
    # ineligible; a future candidate cannot bypass the required provenance,
    # uncertainty, observable, and adversarial-control fields.
    contract = HERE / "data" / "external-comparison.json"
    if not contract.is_file():
        print("EXTERNAL PSI VALIDATION CONTRACT")
        print("INELIGIBLE: no external-comparison manifest with source data is committed")
    # Values printed by the source; they are reference metadata, never fit from
    # the executable's output.
    reference = {
        "nozzle_d_m": 0.0095,
        "stand_off_m": 0.0762,
        "velocity_m_s": 34.0,
        "grain_range_m": (0.0002, 0.0006),
        "duration_s": 100.0,
        "geometry": "circular nitrogen pipe / quartz sand",
    }
    current = {
        "grain_d_m": value("grain_d"),
        "domain_height_m": value("z_hi"),
        "duration_s": value("dyn_steps") * value("sph_dt"),
        "geometry": "planar thin slab / air-like gas",
    }
    failures = [
        "circular nozzle is absent (the current mesh is planar)",
        f"grain diameter {current['grain_d_m']:.4g} m is outside the published "
        f"{reference['grain_range_m'][0]:.4g}–{reference['grain_range_m'][1]:.4g} m range",
        f"run duration {current['duration_s']:.4g} s is not the published {reference['duration_s']:.0f} s depth record",
        "source gas/material (nitrogen/quartz) is not configured",
        "current output counts entrained parcels; it does not measure the published crater-depth trace",
    ]
    print("EXTERNAL PSI REFERENCE AUDIT")
    print("source: Metzger et al. 2009, doi:10.1063/1.3180041, arXiv:0905.4851")
    print("published case:", reference["geometry"], f"d={reference['nozzle_d_m']:.4g} m", f"H={reference['stand_off_m']:.4g} m", f"U={reference['velocity_m_s']:.1f} m/s")
    print("executable case:", current["geometry"], f"grain={current['grain_d_m']:.4g} m", f"duration={current['duration_s']:.4g} s")
    for failure in failures:
        print("INELIGIBLE:", failure)
    print("EXTERNAL PSI COMPARISON: INELIGIBLE (no validation verdict may be emitted)")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
