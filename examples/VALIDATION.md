# Evidence and smoke-execution index

## PSI acceptance status

The pending impinging-plume/crater goal is **not validated by this repository**.
The available smoke case neither advances a CFD jet nor reports a crater,
erosion, or ejecta observable. A future PSI acceptance case must satisfy the
held-out protocol in
[`packed_bed_seam/EXTERNAL_VALIDATION.md`](packed_bed_seam/EXTERNAL_VALIDATION.md).
Until then, the result below is only a packed-bed seam smoke execution.

## Coupled SPH-CFD packed-bed seam smoke case

[The packed_bed_seam smoke case](packed_bed_seam/README.md) runs the live SPH-CFD
**force hand-off** at minimum fluidization, with an imposed homogeneous gas
velocity rather than an advancing CFD calculation. It reports the measured SPH
`U_mf`, a Wen-Yu packed-bed comparator, a resolved same-seam DEM consistency
value, two executable altered-coupling probes (omitted pressure gradient and an
incorrect voidage exponent), and a dynamic sweep above and below `U_mf`.

The check is a finite, coarse smoke execution, not experimental proof of
plume-surface interaction. It deliberately carries no local pass bands:
configuring an expected answer and then accepting the same output is not
independent validation. In particular, the DEM value is same-seam consistency;
Wen--Yu is a packed-bed correlation, not a matched impinging-jet experiment.
