# External validation contract for an impinging-plume claim

`plume_surface` does not meet this contract. Its imposed, homogeneous velocity
field and packed-bed `U_mf` diagnostic are deliberately excluded from plume,
crater, erosion, or ejecta acceptance.

## Admissible evidence

An acceptance submission must include a primary experimental dataset (or its
publisher-approved digitization) with a stable source URL/DOI, acquisition date,
license, raw or traceable measured values, and a cryptographic hash of the exact
input file. The experiment must match the simulated observable and report, at
minimum, nozzle geometry and stand-off, gas composition/state and mass flow or
chamber conditions, gravity/ambient pressure, granular material and preparation,
observation times, and measurement uncertainty.

NASA's PSI program describes controlled plume--surface testing, but its public
program pages and NTRS reports are **not** treated here as numerical data:

- NASA NTRS, *Plume Surface Interaction Physics Focused Ground Test*, 2021,
  https://ntrs.nasa.gov/citations/20210016650
- NASA, *What a Blast: NASA Langley Begins Plume Surface Interaction Tests*,
  https://www.nasa.gov/general/what-a-blast-nasa-langley-begins-plume-surface-interaction-tests/

These sources establish that an appropriate experimental program exists; they do
not supply a matched, machine-readable acceptance series in this repository.
They therefore cannot be used to mark a simulation pass.

## Pre-registered comparison

Before inspecting the held-out observations, freeze the mesh/time-step study,
material parameters, initial packing/preparation, nozzle/boundary conditions,
observable extraction method, and error metric. Split calibration and evaluation
cases by forcing condition or run identifier. Report every held-out case, signed
residuals, uncertainty treatment, convergence evidence, and failures; no
post-hoc tolerance changes or case removal are allowed.

The simulated and observed quantity must be identical (for example crater-depth
versus time at the stated reference plane, or mass-loss rate over the stated
window). A fluidization threshold, visually similar crater, or a different
material/geometry is not a substitute.

## Adversarial control

Run a deliberately wrong but executable coupling variant on the same frozen
cases: at least one of omitted pressure-gradient feedback, disabled particle-to-
gas momentum reaction, or an independently justified wrong drag closure. It
must degrade the pre-registered held-out metric relative to the nominal model.
If it does not, the proposed observable is not discriminating and the acceptance
attempt fails rather than relaxing a threshold.

## Decision boundary

Only a domain reviewer may accept the resulting evidence package. CI may check
file integrity and reproducibility, but a CI pass is never scientific
validation. Until the package above exists and has been reviewed, repository
documentation must state **not validated for plume-surface predictions**.

## Authorship and limits

This protocol and the surrounding code documentation were prepared with AI
assistance. They are a scope and evidence contract, not expert validation or a
substitute for PSI experimental judgment.
