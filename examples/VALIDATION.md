# Evidence and smoke-execution index

## PSI acceptance status

The pending impinging-plume/crater goal is **not validated by this repository**.
The available smoke case neither advances a CFD jet nor reports a crater,
erosion, or ejecta observable. A future PSI acceptance case must satisfy the
held-out protocol in
[`packed_bed_seam/EXTERNAL_VALIDATION.md`](packed_bed_seam/EXTERNAL_VALIDATION.md).
Until then, the result below is only a packed-bed seam smoke execution.

The named external PSI context records have a separate networked citation-
maintenance audit (`python3 references/audit_public_records.py`). It fails
closed if a public record changes or cannot be retrieved, but it is not a
physical acceptance criterion and is intentionally absent from the executable
validation manifest.

## Boundary-driven CFD/SPH jet case (exploratory)

[The jet_crater case](jet_crater/README.md) restores an advancing CFD inflow
and port-coupled SPH surface response.  Its committed eligibility figure and
external-reference audit deliberately report **INELIGIBLE**, not a PSI pass:
the available Metzger et al. crater trace has a different circular geometry,
gas/material, grain range, duration, and observable.  This is executable seam
evidence plus a fail-closed reference check, not acceptance evidence for the
pending plume/crater goal.

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
