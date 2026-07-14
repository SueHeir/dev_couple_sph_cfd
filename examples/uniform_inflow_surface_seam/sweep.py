#!/usr/bin/env python3
"""Run an exploratory seam probe, then fail closed on an ineligible reference.

This driver intentionally has no eligibility manifest or local diagnostic pass
logic. A manifest authored next to a simulation is not independent evidence.
Until a case-matched primary observation series and a quantitative comparator
exist, the only scientifically honest outcome is an ineligibility report.
"""
from __future__ import annotations

import pathlib
import subprocess

import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "plots" / "external_comparison_ineligible.png"


def audit_reference_mismatch() -> str:
    """Run the known-reference mismatch audit; nonzero is expected for now."""
    proc = subprocess.run(
        ["python3", str(HERE / "external_reference_audit.py")],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    print(proc.stdout)
    if proc.returncode == 0:
        raise RuntimeError("an ineligible-reference audit unexpectedly succeeded")
    if "EXTERNAL PSI COMPARISON: INELIGIBLE" not in proc.stdout:
        raise RuntimeError("reference audit did not emit its fail-closed verdict")
    return proc.stdout


def render_ineligible(audit: str) -> None:
    """Make the non-comparability visible without manufacturing a pass band."""
    reasons = [line.removeprefix("INELIGIBLE: ") for line in audit.splitlines()
               if line.startswith("INELIGIBLE:")]
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    ax.set_axis_off()
    ax.text(0.5, 0.88, "External PSI comparison: INELIGIBLE", ha="center", va="center",
            color="#991b1b", fontsize=20, fontweight="bold")
    ax.text(0.5, 0.76,
            "Metzger et al. (2009) provides crater-depth data, but this executable\n"
            "may not compare its parcel-count diagnostic to that trace.",
            ha="center", va="center", fontsize=11)
    ax.text(0.07, 0.62, "Reference/case mismatches", fontweight="bold", fontsize=13)
    for i, reason in enumerate(reasons):
        ax.text(0.09, 0.54 - i * 0.075, f"• {reason}", fontsize=10, va="top")
    ax.text(0.5, 0.10,
            "No PSI-validation verdict is emitted. A future comparison needs matched geometry,\n"
            "materials, forcing and duration; a same-observable source series with uncertainty;\n"
            "convergence evidence; and an adversarial control against that held-out series.",
            ha="center", va="center", fontsize=10, color="#374151")
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=160)


def main() -> None:
    # Exercise the boundary-driven CFD→SPH path before the reference audit. The
    # audit must never hide a broken executable path, but it also must never
    # turn this path into a scientific acceptance route.
    proc = subprocess.run(
        ["cargo", "run", "--release", "--example", "uniform_inflow_surface_seam", "--",
         "examples/uniform_inflow_surface_seam/config.toml"],
        cwd=ROOT, text=True, check=False,
    )
    if proc.returncode:
        raise SystemExit(proc.returncode)
    audit = audit_reference_mismatch()
    render_ineligible(audit)
    print(f"Wrote {OUT.relative_to(ROOT)} (external PSI comparison ineligible)")
    raise SystemExit("VALIDATION: FAIL — seam executed, but PSI comparison is ineligible")


if __name__ == "__main__":
    main()
