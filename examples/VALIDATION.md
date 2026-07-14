# Validation index

## PSI acceptance status

The pending impinging-plume/crater goal is **not validated by this repository**.
The available regression neither advances a CFD jet nor reports a crater,
erosion, or ejecta observable. A future PSI acceptance case must supply a
geometry-, material-, gas-, forcing-duration-, and observation-time-matched
primary-source series; an advancing CFD case with declared boundaries; and a
same-observable comparator whose deliberately wrong coupling fails. Until then,
the result below is only a packed-bed seam regression.

## Coupled SPH-CFD plume/surface minimum fluidization

[The plume_surface regression](plume_surface/README.md) runs the live SPH-CFD
**force hand-off** at minimum fluidization, with an imposed homogeneous gas
velocity rather than an advancing CFD calculation. Its committed
[result figure](plume_surface/plots/plume_surface_validation.svg) compares the
measured SPH `U_mf` with the independent Wen-Yu correlation and a resolved
DEM-CFD cross-method reference, with their tolerance bands shown. It also
shows two executable fault controls (omitted pressure gradient and an
incorrect voidage exponent) outside those bands, plus the dynamic
pressure/onset gate above and below `U_mf`.

The check is a finite, coarse seam regression, not experimental proof
of plume-surface interaction. In particular, the DEM value is same-seam
consistency; Wen--Yu is a packed-bed correlation, not a matched impinging-jet
experiment. Its reference and model limitations are stated with the case
configuration and source output.
