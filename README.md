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

## What it does — plume ↔ surface interaction

Reusable SPH-CFD coupling code lives in `crates/sph_cfd`: packed-bed closure and
reference helpers, parcel deposition, force-balance measurement, SPH force import,
CFD-side seam resources/systems, and the standard `grass_multi` exchange schedule.
Examples keep case geometry, validation tolerances, comparison packings, and plots.

The `plume_surface` example couples a compressible gas jet (dev_field_efvm, on a
FIELD mesh) to a granular bed (dev_soil_sph, as SOIL particles): the gas exerts
drag on the grains and the grains displace/block the gas — the dynamic
minimum-fluidization / landing-plume-on-regolith cratering problem. The two
solvers run as **grass sub-Apps under one parent schedule** (`Tick → Couple`),
sharing exactly one `grass_app::App` and `soil_core::Atom` type across the seam.

```bash
# all partner repos are sibling checkouts (grass, soil, field, dev_soil_sph, dev_field_efvm)
cargo run --release --example plume_surface -- examples/plume_surface/config.toml
```

Repo-level validation is declared in `validation/manifest.toml` and runs through
`validation/run.sh`, which delegates to the plume-surface sweep so the numeric
gate and committed figure come from the same measured output.

## License

MIT OR Apache-2.0
