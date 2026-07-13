# PSI external-evidence gate

This repository has a minimum-fluidization seam check in
[`../plume_surface/`](../plume_surface/). That check is not an impinging-plume
crater, erosion-rate, or ejecta validation.

The associated goal requires a primary-source data series that matches nozzle
geometry, stand-off, gas, grain material/range, forcing duration, and the
reported crater observable. Metzger et al., *Craters Formed in Granular Beds by
Impinging Jets of Gas* (2009),
[doi:10.1063/1.3180041](https://doi.org/10.1063/1.3180041), is relevant context,
but it is not a usable reference here: this repository has no digitized,
condition-matched data or corresponding advancing-CFD case.

`audit.py` intentionally exits nonzero until a manifest supplies a cited source
table with checksum and measurement uncertainty, independently generated model
predictions, a separately generated wrong-coupling control, and a comparator.
It cannot turn an unmatched publication, an analytic inlet profile, or the
fluidization figure into PSI-validation evidence.

```bash
python3 examples/psi_evidence/audit.py
```

AI-assisted authorship disclosure: this gate and its documentation were prepared
with AI assistance. No new experiment, digitization, uncertainty estimate, or
PSI validation result is claimed.
