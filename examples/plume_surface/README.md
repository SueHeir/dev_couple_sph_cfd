# plume_surface

This example checks an SPH-to-CFD **packed-bed seam regression** at the dynamic
minimum-fluidization limit. Its gas velocity is imposed every coupling tick;
`CfdStatePlugin` and `IdealGasPlugin` are present as the data carrier, but no
`SolverPlugin` advances a CFD field. It measures the SPH-continuum `U_mf` from
the live force hand-off, compares it with Wen-Yu, runs two fault controls that
must move outside tolerance, and sweeps the coupled bed through onset.

![plume_surface validation result](plots/plume_surface_validation.svg)

The plot is generated from `sweep.py`, which runs the example and parses its own
reported `U_mf` values, tolerance checks, negative controls, and dynamic pressure
sweep. The current committed figure shows `VALIDATION: PASS` for this limited
seam regression only. Its DEM comparison is same-seam consistency, not
independent PSI evidence. It is not evidence for an advancing impinging-jet
crater, erosion-rate, or ejecta prediction; see
[`../psi_evidence/`](../psi_evidence/) for the fail-closed external-evidence
requirements for that separate claim.

```bash
python3 examples/plume_surface/sweep.py
```
