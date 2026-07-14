# Independent reference provenance

- L. Schiller and A. Naumann, *Z. Ver. Deut. Ing.* 77, 318 (1935), sphere-drag closure used for the local coupling force.

This reference supplies a local drag closure only. It is not a plume-surface
experiment and is not used as an external PSI target. The severed-port run is a
seam observation, not a substitute for independent data.

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
