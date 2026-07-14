# Evidence and regression index

## PSI acceptance status

The pending impinging-plume/crater goal is **not validated by this repository**.
The available regression neither advances a CFD jet nor reports a crater,
erosion, or ejecta observable. A future PSI acceptance case must satisfy the
held-out protocol in
[`packed_bed_seam/EXTERNAL_VALIDATION.md`](packed_bed_seam/EXTERNAL_VALIDATION.md).
Until then, the result below is only a packed-bed seam regression.

## Coupled SPH-CFD packed-bed seam at minimum fluidization

[The packed_bed_seam regression](packed_bed_seam/README.md) runs the live SPH-CFD
**force hand-off** at minimum fluidization, with an imposed homogeneous gas
velocity rather than an advancing CFD calculation. Its committed
[result figure](packed_bed_seam/plots/packed_bed_seam_validation.svg) compares the
measured SPH `U_mf` with a Wen-Yu packed-bed comparator and a resolved
same-seam DEM consistency value, with frozen regression bands shown. It also
shows two executable fault controls (omitted pressure gradient and an
incorrect voidage exponent) outside those bands, plus the dynamic
pressure/onset gate above and below `U_mf`.

The check is a finite, coarse seam regression, not experimental proof
of plume-surface interaction. The bands are regression limits, not confidence
intervals or experimental tolerances. In particular, the DEM value is same-seam
consistency; Wen--Yu is a packed-bed correlation, not a matched impinging-jet
experiment.
