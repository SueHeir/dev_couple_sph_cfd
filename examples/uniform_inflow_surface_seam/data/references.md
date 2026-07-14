# Independent reference provenance

- R. A. Bagnold, *The Physics of Blown Sand and Desert Dunes* (1941), aerodynamic entrainment scaling `u* = A sqrt(((rho_s-rho_f)/rho_f) g d)`, with `A` of order 0.1.
- J. D. Iversen and B. R. White, “Saltation threshold on Earth, Mars and Venus,” *Sedimentology* 29, 111–119 (1982), experimental/aerodynamic threshold context used for the deliberately broad `A=[0.06,0.20]` acceptance band.
- L. Roberts, “The action of a hypersonic jet on a dust layer,” IAS Paper 63-50 (1963), wall-jet plume/surface erosion and off-axis shear mechanism.
- L. Schiller and A. Naumann, *Z. Ver. Deut. Ing.* 77, 318 (1935), independent sphere-drag closure used to recover onset; it does not contain the Bagnold coefficient.

These references motivate closures and qualitative mechanisms only. This case does not include a digitized Roberts crater/erosion series at matched nozzle, stand-off, gas, grain, and gravity conditions, so neither the Bagnold band nor off-axis response is an external PSI validation gate. The cohesive and severed-port runs are implementation controls, not substitutes for independent data.

Metzger et al. (2009), doi:10.1063/1.3180041 (open preprint arXiv:0905.4851),
is now retained as an independently checked crater-depth source.  Its case is
not comparable to this executable; see [metzger_2009_reference.md](metzger_2009_reference.md).
The executable `external_reference_audit.py` fails closed on that mismatch.

The example was drafted with AI assistance. This rescue revision replaces the
analytic interior wall-jet initialization with a quiescent gas field and a physical
downward `ZHi` inflow boundary, which the CFD solver advances to the bed. That is a
software/coupling correction, not scientific validation. Its current limit is
explicit: it does not establish quantitative plume-surface accuracy. A future
validation must commit source values and units, compare a measured
crater-depth/volume/ejecta-rate observable with a predeclared defensible tolerance,
and include an adversarial control that fails that external comparison.

Those minimum evidence fields are review requirements, not a local manifest
that can authorize a verdict. A future comparator must independently reproduce
a case-matched primary observation series and evaluate both the model and a
wrong-coupling control against it.

AI-assisted implementation and the absence of an independently reproduced,
geometry-matched experimental series remain validation limitations.
