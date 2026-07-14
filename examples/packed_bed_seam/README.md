# packed_bed_seam

This example executes an SPH-to-CFD **packed-bed seam smoke case** at the dynamic
minimum-fluidization limit. Its gas velocity is imposed every coupling tick;
`CfdStatePlugin` and `IdealGasPlugin` are present as the data carrier, but no
`SolverPlugin` advances a CFD field. It reports the SPH-continuum `U_mf` from
the live force hand-off, a same-seam DEM value, a Wen-Yu comparator, two
altered-coupling sensitivity probes, and a dynamic sweep. It deliberately
applies no local numerical acceptance thresholds: these values are diagnostics,
not an independent test. Successful completion means only that the configured
local execution path ran. The DEM comparison shares the seam, while Wen--Yu is
a packed-bed correlation; neither is independent PSI evidence. It is not
evidence for an advancing impinging-jet
crater, erosion-rate, or ejecta prediction. The only admissible route to that
claim is the held-out, adversarial protocol in
[`EXTERNAL_VALIDATION.md`](EXTERNAL_VALIDATION.md).

```bash
cargo run --release --example packed_bed_seam -- examples/packed_bed_seam/config.toml
```
