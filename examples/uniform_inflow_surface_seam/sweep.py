#!/usr/bin/env python3
"""Render exploratory seam controls, but refuse to label them PSI validation.

The executable case is deliberately audited before any plot is written.  While
the available external crater trace is geometry-ineligible, this driver creates
an eligibility figure instead of a deceptively green validation chart.
"""
from __future__ import annotations

import pathlib
import re
import subprocess

import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "plots" / "external_comparison_ineligible.png"


def external_comparison_is_eligible() -> tuple[bool, str]:
    """Run the independent-reference firewall; nonzero is expected for now."""
    contract = subprocess.run(
        ["python3", str(HERE / "data" / "validation_contract.py"),
        str(HERE / "data" / "external-comparison.json")],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    print(contract.stdout)
    proc = subprocess.run(
        ["python3", str(HERE / "external_reference_audit.py")],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    print(proc.stdout)
    audit = contract.stdout + proc.stdout
    eligible = (contract.returncode == 0 and proc.returncode == 0
                and "EXTERNAL PSI COMPARISON: ELIGIBLE" in proc.stdout)
    # A zero exit without the explicit eligibility marker is not authorization to
    # compare results.  Conversely, the known mismatch must be fail-closed.
    if proc.returncode == 0 and contract.returncode == 0 and not eligible:
        raise RuntimeError("reference audit returned success without eligibility marker")
    return eligible, audit


def render_ineligible(audit: str) -> None:
    """Commit a result figure that makes the missing external comparison visible."""
    reasons = [line.removeprefix("INELIGIBLE: ") for line in audit.splitlines()
               if line.startswith("INELIGIBLE:")]
    fig, ax = plt.subplots(figsize=(10, 5.2), constrained_layout=True)
    ax.set_axis_off()
    ax.text(0.5, 0.88, "External PSI comparison: INELIGIBLE", ha="center", va="center",
            color="#991b1b", fontsize=20, fontweight="bold")
    ax.text(0.5, 0.76,
            "Metzger et al. (2009), doi:10.1063/1.3180041, provides crater-depth data;\n"
            "the executable must not compare its parcel-count diagnostic to that trace.",
            ha="center", va="center", fontsize=11)
    ax.text(0.07, 0.62, "Independent audit failures", fontweight="bold", fontsize=13)
    for i, reason in enumerate(reasons):
        ax.text(0.09, 0.54 - i * 0.075, f"• {reason}", fontsize=10, va="top")
    ax.text(0.5, 0.10,
            "No PASS/FAIL PSI-validation verdict is emitted.  The live CFD→SPH controls remain exploratory.\n"
            "A future comparison needs matched circular geometry, materials, duration, a crater-depth observable,\n"
            "source data with uncertainty, and an adversarial external-comparison control.",
            ha="center", va="center", fontsize=10, color="#374151")
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=160)


def grab(pattern: str, text: str) -> re.Match[str]:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"missing validation datum: {pattern}")
    return match


def main() -> None:
    # Always exercise the boundary-driven CFD→SPH path.  The reference audit below
    # decides only whether this execution may be compared to a PSI observation;
    # it must not become a shortcut that hides a broken runnable case.
    proc = subprocess.run(
        ["cargo", "run", "--release", "--example", "uniform_inflow_surface_seam", "--", "examples/uniform_inflow_surface_seam/config.toml"],
        cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    print(proc.stdout)
    if proc.returncode:
        raise SystemExit(proc.returncode)
    a = grab(r"Bagnold A_meas=([0-9.]+).*band \[([0-9.]+),([0-9.]+)\]", proc.stdout)
    exponent = grab(r"exponent p=([0-9.]+).*cohesive-control p=([0-9.]+)", proc.stdout)
    rows = re.findall(r"^\s*([0-9.]+)\s+([0-9.]+)\s+(\d+)/(\d+).*?\s+([0-9.]+)\s+(?:packed|ERODES)", proc.stdout, re.MULTILINE)
    if not rows:
        raise RuntimeError("missing live erosion sweep")
    factors = [float(row[0]) for row in rows]
    eroded = [int(row[2]) for row in rows]
    surface = [int(row[3]) for row in rows]
    offsets = [float(row[4]) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    ax = axes[0]
    ref, lo, hi = 0.10, float(a.group(2)), float(a.group(3))
    ax.axhspan(lo, hi, color="#dbeafe", label="Bagnold/Iversen–White band")
    ax.axhline(ref, color="#7c3aed", ls="--", label="published A≈0.10")
    ax.scatter(["coupled onset"], [float(a.group(1))], color="#2563eb", s=65, zorder=3)
    ax.set_ylim(0, max(0.23, hi * 1.15)); ax.set_ylabel("Bagnold coefficient A")
    ax.set_title(f"Mechanism diagnostic: {verdict}\n g exponent {exponent.group(1)}; cohesive control {exponent.group(2)}")
    ax.legend(fontsize=8, loc="upper right")
    ax = axes[1]
    ax.plot(factors, eroded, "o-", color="#dc2626", label="entrained surface parcels")
    ax.plot(factors, surface, "o--", color="#64748b", label="surface parcels")
    for x, y, off in zip(factors, eroded, offsets): ax.annotate(f"offset/a={off:.2f}", (x, y), xytext=(3, 6), textcoords="offset points", fontsize=8)
    ax.axvline(1, color="#7c3aed", ls="--", label="onset")
    ax.set(xlabel="uniform-inflow strength U/u_gc", ylabel="parcel count", title="Live CFD→SPH exploratory response")
    ax.legend(fontsize=8)
    fig.suptitle("Uniform-inflow seam diagnostics — not external PSI validation", fontweight="bold")

    eligible, audit = external_comparison_is_eligible()
    if not eligible:
        # Do not overwrite the committed eligibility figure with local diagnostic
        # curves: a reader must see why these curves cannot be called a PSI match.
        plt.close(fig)
        render_ineligible(audit)
        print(f"Wrote {OUT.relative_to(ROOT)} (external PSI comparison ineligible)")
        raise SystemExit("VALIDATION: FAIL — CFD/SPH seam executed, but external PSI comparison is ineligible")

    # An eligible source still needs a comparator with source uncertainty and a
    # wrong-coupling run.  Never convert eligibility alone into success.
    plt.close(fig)
    raise RuntimeError("eligible reference has no implemented quantitative comparator")


if __name__ == "__main__":
    main()
