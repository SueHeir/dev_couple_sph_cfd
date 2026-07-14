#!/usr/bin/env python3
"""Run the packed-bed regression and render diagnostic output."""

from __future__ import annotations

import html
import pathlib
import re
import subprocess
import sys


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONFIG = HERE / "config.toml"
PLOTS = HERE / "plots"
FIGURE = PLOTS / "plume_surface_validation.svg"


def run_case() -> str:
    cmd = [
        "cargo",
        "run",
        "--release",
        "--example",
        "plume_surface",
        "--",
        str(CONFIG.relative_to(ROOT)),
    ]
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    sys.stdout.write(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(proc.returncode)
    return proc.stdout


def must(pattern: str, text: str, label: str) -> re.Match[str]:
    m = re.search(pattern, text, re.MULTILINE)
    if not m:
        raise SystemExit(f"could not parse {label}")
    return m


def parse_output(text: str) -> dict[str, object]:
    sph = float(must(r"SPH continuum .* : ([0-9.]+)", text, "SPH U_mf").group(1))
    wen_m = must(
        r"Wen & Yu 1966 correlation .* : ([0-9.]+)\s+rel.err ([0-9.]+)%\s+\(regression limit ([0-9.]+)%\)",
        text,
        "Wen-Yu U_mf",
    )
    dem_m = must(
        r"DEM discrete .* : ([0-9.]+)\s+consistency rel.err ([0-9.]+)%\s+\(regression limit ([0-9.]+)%\)",
        text,
        "DEM-CFD U_mf",
    )
    neg_m = must(
        r"fault controls: omit-gradP ([0-9.]+) \(([+-][0-9.]+)%\)\s+eps-power-bug ([0-9.]+) \(([+-][0-9.]+)%\)",
        text,
        "fault controls",
    )
    result = must(r"REGRESSION: (PASS|FAIL).*", text, "regression verdict").group(0)

    sweep = []
    row_re = re.compile(
        r"^\s*([0-9.]+)\s+([0-9.]+)\s+([+-][0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.]+)\s+(.*)$",
        re.MULTILINE,
    )
    for row in row_re.finditer(text):
        sweep.append(
            {
                "factor": float(row.group(1)),
                "u": float(row.group(2)),
                "vz": float(row.group(3)),
                "pressure": float(row.group(4)),
                "pfrac": float(row.group(5)),
                "state": row.group(6).strip(),
            }
        )
    if not sweep:
        raise SystemExit("could not parse dynamic pressure sweep")

    return {
        "sph": sph,
        "wenyu": float(wen_m.group(1)),
        "wenyu_err": float(wen_m.group(2)),
        "tol_wenyu": float(wen_m.group(3)),
        "dem": float(dem_m.group(1)),
        "dem_err": float(dem_m.group(2)),
        "tol_dem": float(dem_m.group(3)),
        "omit_gradp": float(neg_m.group(1)),
        "omit_gradp_shift": float(neg_m.group(2)),
        "eps_bug": float(neg_m.group(3)),
        "eps_bug_shift": float(neg_m.group(4)),
        "sweep": sweep,
        "result": result,
    }


def sx(x: float, xmin: float, xmax: float, left: float, width: float) -> float:
    return left + (x - xmin) / (xmax - xmin) * width


def sy(y: float, ymin: float, ymax: float, top: float, height: float) -> float:
    return top + height - (y - ymin) / (ymax - ymin) * height


def line(x1: float, y1: float, x2: float, y2: float, color: str, width: float = 2.0, dash: str = "") -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{width}"{dash_attr}/>'


def text(x: float, y: float, s: str, size: int = 14, anchor: str = "start", weight: str = "400") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="#1f2933">{html.escape(s)}</text>'
    )


def render_svg(data: dict[str, object]) -> str:
    w, h = 1120, 640
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif}.small{font-size:12px;fill:#52616f}</style>',
        text(36, 34, "packed-bed force-transfer regression: diagnostics and fault controls", 20, weight="700"),
        text(36, 58, str(data["result"]), 13),
    ]

    left, top, width, height = 74, 98, 430, 410
    sph = float(data["sph"])
    dem = float(data["dem"])
    wenyu = float(data["wenyu"])
    omit_gradp = float(data["omit_gradp"])
    eps_bug = float(data["eps_bug"])
    tol_dem = float(data["tol_dem"]) / 100.0
    tol_wenyu = float(data["tol_wenyu"]) / 100.0
    vals = [sph, dem, wenyu, omit_gradp, eps_bug]
    xmin, xmax = 0.0, max(vals) * 1.12

    parts += [
        text(left, top - 20, "Measured U_mf vs references and negative controls", 15, weight="700"),
        line(left, top, left, top + height, "#9aa5b1", 1),
        line(left, top + height, left + width, top + height, "#9aa5b1", 1),
    ]
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = left + frac * width
        v = xmin + frac * (xmax - xmin)
        parts.append(line(x, top, x, top + height, "#e4e7eb", 1))
        parts.append(text(x, top + height + 24, f"{v:.3f}", 11, "middle"))
    parts.append(text(left + width / 2, top + height + 52, "superficial velocity U_mf [m/s]", 13, "middle"))

    # Regression bands are frozen local guardrails, not uncertainty or validation bounds.
    for ref, tol, color, y0, label in [
        (dem, tol_dem, "#d8f3dc", top + 54, f"same-seam regression +/- {data['tol_dem']:.0f}%"),
        (wenyu, tol_wenyu, "#dbeafe", top + 138, f"Wen-Yu regression +/- {data['tol_wenyu']:.0f}%"),
    ]:
        x0 = sx(ref * (1.0 - tol), xmin, xmax, left, width)
        x1 = sx(ref * (1.0 + tol), xmin, xmax, left, width)
        parts.append(f'<rect x="{x0:.1f}" y="{y0 - 26:.1f}" width="{x1 - x0:.1f}" height="52" fill="{color}" opacity="0.9"/>')
        parts.append(text(x1 + 6, y0 + 4, label, 11))

    rows = [
        ("SPH measured", sph, "#2f80ed"),
        ("same-seam DEM", dem, "#219653"),
        ("Wen-Yu comparator", wenyu, "#8e44ad"),
        ("omit grad-P", omit_gradp, "#d64545"),
        ("eps-power bug", eps_bug, "#d64545"),
    ]
    for i, (label, value, color) in enumerate(rows):
        y = top + 42 + i * 78
        x = sx(value, xmin, xmax, left, width)
        parts.append(line(left, y, x, y, color, 13))
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="{color}"/>')
        parts.append(text(left - 10, y + 5, label, 12, "end"))
        parts.append(text(x + 10, y + 5, f"{value:.4f}", 12))

    right, rtop, rwidth, rheight = 620, 98, 430, 410
    sweep = list(data["sweep"])  # type: ignore[arg-type]
    xs = [float(row["factor"]) for row in sweep]
    ys = [float(row["pfrac"]) for row in sweep]
    xmin2, xmax2 = min(xs) * 0.9, max(xs) * 1.05
    ymin2, ymax2 = 0.0, max(ys + [0.12]) * 1.12
    parts += [
        text(right, rtop - 20, "Dynamic pressure sweep through onset", 15, weight="700"),
        line(right, rtop, right, rtop + rheight, "#9aa5b1", 1),
        line(right, rtop + rheight, right + rwidth, rtop + rheight, "#9aa5b1", 1),
    ]
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        yv = ymin2 + frac * (ymax2 - ymin2)
        y = sy(yv, ymin2, ymax2, rtop, rheight)
        parts.append(line(right, y, right + rwidth, y, "#e4e7eb", 1))
        parts.append(text(right - 10, y + 4, f"{yv:.2f}", 11, "end"))
    for xv in [0.5, 0.8, 1.0, 1.2, 1.5]:
        if xmin2 <= xv <= xmax2:
            x = sx(xv, xmin2, xmax2, right, rwidth)
            parts.append(line(x, rtop, x, rtop + rheight, "#eef2f7", 1))
            parts.append(text(x, rtop + rheight + 24, f"{xv:.1f}", 11, "middle"))
    onset_x = sx(1.0, xmin2, xmax2, right, rwidth)
    fluid_y = sy(0.05, ymin2, ymax2, rtop, rheight)
    packed_y = sy(0.10, ymin2, ymax2, rtop, rheight)
    parts += [
        line(onset_x, rtop, onset_x, rtop + rheight, "#ef8354", 2, "6 5"),
        text(onset_x + 6, rtop + 18, "U/U_mf = 1", 11),
        line(right, fluid_y, right + rwidth, fluid_y, "#2f80ed", 2, "6 5"),
        text(right + rwidth - 4, fluid_y - 7, "fluidized p/p0 <= 0.05", 11, "end"),
        line(right, packed_y, right + rwidth, packed_y, "#219653", 2, "6 5"),
        text(right + rwidth - 4, packed_y + 16, "packed p/p0 > 0.10", 11, "end"),
    ]
    pts = [
        (sx(float(row["factor"]), xmin2, xmax2, right, rwidth), sy(float(row["pfrac"]), ymin2, ymax2, rtop, rheight))
        for row in sweep
    ]
    if len(pts) > 1:
        d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
        parts.append(f'<polyline points="{d}" fill="none" stroke="#111827" stroke-width="3"/>')
    for row, (x, y) in zip(sweep, pts):
        color = "#219653" if float(row["factor"]) < 1.0 else "#2f80ed"
        parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
        parts.append(text(x + 8, y - 8, f"{float(row['pfrac']):.3f}", 11))
    parts.append(text(right + rwidth / 2, rtop + rheight + 52, "superficial velocity factor U/U_mf", 13, "middle"))
    parts.append(text(right - 48, rtop + rheight / 2, "contact pressure fraction p/p0", 13, "middle"))

    parts.append(text(36, 602, f"Fault-control shifts vs Wen-Yu: omit grad-P {data['omit_gradp_shift']:+.1f}%, eps-power bug {data['eps_bug_shift']:+.1f}%.", 13))
    parts.append("</svg>")
    return "\n".join(parts)


def main() -> None:
    output = run_case()
    data = parse_output(output)
    PLOTS.mkdir(exist_ok=True)
    FIGURE.write_text(render_svg(data), encoding="utf-8")
    print(f"Wrote {FIGURE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
