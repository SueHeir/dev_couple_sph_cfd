# plume_surface

This example validates the coupled SPH-CFD plume/surface seam at the dynamic
minimum-fluidization limit. It measures the SPH-continuum `U_mf` from the live
coupling, compares it with the same-seam DEM-CFD reference and the Wen-Yu
correlation, runs two negative controls that must move outside tolerance, and
sweeps the live coupled bed through onset.

![plume_surface validation result](plots/plume_surface_validation.svg)

The plot is generated from `sweep.py`, which runs the example and parses its own
reported `U_mf` values, tolerance checks, negative controls, and dynamic pressure
sweep. The current committed figure shows `VALIDATION: PASS`.

```bash
python3 examples/plume_surface/sweep.py
```
