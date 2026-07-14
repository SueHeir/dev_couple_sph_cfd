# dev_couple_sph_cfd

<!-- disclaimer-banner -->
> This code was fully written via **Claude 4.6,4.8 and Fable 5**, and stands as a proof of concept for a **bevy-like** ecosystem for physics simulation research, with the goal of testing if one scheduler/framework (**GRASS**) works for most scientific codes. **SOIL** and **FIELD** are particle- and mesh-based substrates for physics such as **DIRT** (DEM) or **dev_field_efvm**. Note that all other physics based repos I have start with **dev_**, as I do **NOT** know these methods. Please read, evaluate, use with a grain of salt, I have not personally read or reviewed everything here.
<!-- /disclaimer-banner -->


A **cross-substrate coupling**: it joins the granular-SPH tier
[dev_soil_sph](https://github.com/SueHeir/dev_soil_sph) (a particle/SOIL solver)
to the compressible-CFD tier
[dev_field_efvm](https://github.com/SueHeir/dev_field_efvm) (a mesh/FIELD solver)
through **GRASS's open-box coupling layer** (`grass_multi`). It is not a physics
tier of its own — it owns no new solver, only the seam between two.

```
GRASS   framework: App, Plugin, Scheduler, coupling (grass_multi)
  ├─ SOIL  (particle substrate) ── dev_soil_sph   granular SPH  (μ(I) continuum)  ┐
  └─ FIELD (mesh substrate)     ── dev_field_efvm  compressible CFD (Riemann/IBM) ┘
                                          └── dev_couple_sph_cfd  ← the coupling (you are here)
```

## Why a separate repo

A coupling that depends on **two** substrate tiers does not belong inside either
one — burying it in `dev_field_efvm` made that CFD tier drag in an SPH dependency
it otherwise has no business with. Cross-substrate couplings (SOIL ↔ FIELD) are
their own thing: they need `grass_multi` + a transport/interphase seam, and they
compose two independently-developed tiers. So each such coupling gets its own
`dev_couple_*` repo, depending on its two partner tiers and nothing more.

## What it does — packed-bed fluidization seam

Reusable SPH-CFD coupling code lives in `crates/sph_cfd`: packed-bed closure and
reference helpers, parcel deposition, force-balance measurement, SPH force import,
CFD-side seam resources/systems, and the standard `grass_multi` exchange schedule.
Examples keep case geometry and diagnostic comparison packings.

The `packed_bed_seam` example couples an **imposed homogeneous interstitial gas
velocity** (stored in dev_field_efvm's FIELD state) to a granular bed
(dev_soil_sph, as SOIL particles). It exercises the drag/pressure-gradient
force hand-off and runs an executable packed-bed seam smoke case with a Wen--Yu
diagnostic comparator and sensitivity probes. It does not advance a CFD
solution or represent a nozzle, a plume, or a crater. The two sub-Apps run
as **grass sub-Apps under one parent schedule** (`Tick → Couple`), sharing
exactly one `grass_app::App` and `soil_core::Atom` type across the seam.

This repository does not currently provide an advancing-CFD impinging-jet case,
nor a validated crater, erosion, or ejecta prediction. Such a claim requires an
advancing CFD solver with specified inlet/outlet conditions, a geometry-,
material-, forcing-, and observation-time-matched experimental series, and a
same-observable comparator with an adversarial wrong-coupling control. A
self-consistent flow profile or differently configured experiment is not a
substitute. The concrete held-out protocol, including source provenance and an
adversarial wrong-coupling control, is in
[`examples/packed_bed_seam/EXTERNAL_VALIDATION.md`](examples/packed_bed_seam/EXTERNAL_VALIDATION.md).
No local manifest or claim guard is scientific evidence. This code and its
documentation were authored with AI assistance and require domain-expert review
before use in scientific conclusions.

The coupling system runs on the **parent** App. It obtains stable participant
handles from `SubApps`, reads each solver's resources between child ticks, and
returns the force through the public seam resources in
[`crates/sph_cfd/src/seam.rs`](crates/sph_cfd/src/seam.rs). Child systems use
ordinary `Res`/`ResMut` for their own state; no child contains a `MultiRes` view
of another solver.

This repo currently proves the local composition path. It does **not** yet claim
the full plume solver under split MPI, but it now proves the same distributed
coupling contract with the runnable
[`routed_sph_cfd`](examples/routed_sph_cfd/main.rs) example: one binary, one
`mode = "auto"` TOML, coupling-owned parcel/force records, and FIELD-owned
position-to-rank lookup over GRASS's generic addressed exchange.

```bash
cargo run --example routed_sph_cfd --features mpi-routing
cargo build --example routed_sph_cfd --features mpi-routing
mpirun --oversubscribe -np 5 target/debug/examples/routed_sph_cfd
```

```bash
# all partner repos are sibling checkouts (grass, soil, field, dev_soil_sph, dev_field_efvm)
cargo run --release --example packed_bed_seam -- examples/packed_bed_seam/config.toml
```

The runnable packed-bed smoke case is intentionally kept separate from the
unmet impinging-plume acceptance claim. It demonstrates a configured local
force hand-off, but it is not an eligibility gate for a crater or erosion
prediction.

## License

MIT OR Apache-2.0
