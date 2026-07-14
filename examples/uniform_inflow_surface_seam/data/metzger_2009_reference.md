# Metzger et al. (2009) external-reference audit

**Source.** P. T. Metzger, R. C. Latta III, J. M. Schuler, and C. D. Immer,
“Craters Formed in Granular Beds by Impinging Jets of Gas,” *AIP Conference
Proceedings* **1145**, 767–772 (2009), doi:10.1063/1.3180041,
arXiv:0905.4851.  This open preprint was retrieved on 2026-07-11.  The
information below is transcribed from its Fig. 2/Fig. 3 caption and the
paragraph immediately preceding Fig. 2, not fitted to this implementation.

## Published case available for comparison

The paper reports a **circular nitrogen** jet from a straight pipe, with a
0.95 cm pipe diameter, 7.62 cm stand-off, and 34 m/s jet velocity, impinging
on sieved quartz construction sand (200–280, 280–300, 300–450, and 500–600
µm).  Its Fig. 3 measures overall and inner crater depth through 100 s.  It
reports logarithmic-period depth growth, `D(t) = a ln(b t)`, for those
experiments.

## Why this is not an eligible reference for the current executable case

`uniform_inflow_surface_seam/config.toml` describes an air-like, planar thin-slab calculation,
with 3 mm grains and 0.016 s of coupled evolution.  It has neither the
circular nozzle nor the source material/diameter/stand-off/time regime.  A
numerical comparison to the Fig. 3 depth trace would therefore be a
cross-geometry extrapolation, not a validation.  The companion
`external_reference_audit.py` is deliberately fail-closed and reports this
mismatch as an **ineligible** comparison.  It must not be inverted into a
passing criterion.

## Required next evidence

Before an external PSI claim can be made, commit a source table (or a
documented digitization with uncertainty) for one experiment and configure
the executable with its nozzle geometry, gas, granular material, grain range,
stand-off, forcing history, and observation interval.  The output must then
measure the same crater-depth/volume/ejecta observable and compare it with a
predeclared uncertainty-supported tolerance, alongside a deliberately wrong
coupling control that fails that same external comparison.
