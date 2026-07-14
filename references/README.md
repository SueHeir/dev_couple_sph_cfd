# Public-reference maintenance

`audit_public_records.py` checks that the two public PSI context pages cited by
the held-out protocol still identify themselves as the expected records. It is
useful for catching a moved or replaced citation.

It does not retrieve experimental observations, compare a simulated observable,
set a tolerance, or make a scientific verdict. It must remain outside
`validation/manifest.toml`; success is not evidence for plume-surface
prediction. The admissible evidence and adversarial-control requirements are
in [`../examples/packed_bed_seam/EXTERNAL_VALIDATION.md`](../examples/packed_bed_seam/EXTERNAL_VALIDATION.md).
