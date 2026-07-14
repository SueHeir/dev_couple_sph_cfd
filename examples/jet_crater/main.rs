//! **Jet–surface interaction: an advancing CFD impinging inflow driving a granular-
//! SPH free surface through the GRASS coupling seam, with exploratory mechanism
//! diagnostics.**
//!
//! A rocket exhaust plume impinging on a granular surface can turn into a wall jet
//! whose near-surface shear entrains surface grains and excavates a **crater** — the
//! plume-surface interaction (PSI) problem central to lunar/Mars landing. This
//! example starts the gas quiescent and supplies a downward top-boundary inflow;
//! the **CFD gas solver** (test-cfd FIELD substrate,
//! `CfdState`) to the **granular-SPH μ(I) continuum** (`sph_core`) *through the
//! grass coupling contract* — `grass_multi`'s exchange **ports** (`add_port`,
//! `expose_field`, `consume_field`) — and advances the coupled erosion state. It
//! is deliberately *not* a quantitative PSI validation: this repository has no
//! digitized, geometry-matched crater or erosion-rate data. Its controls test the
//! executable seam, not agreement with a published crater data series.
//!
//! ## What is coupled, and through what seam (the goal)
//!
//! The two solvers share **only the port contract types** — neither names the
//! other's namespace or internal resources:
//!   - the SPH bed **exposes** its free-surface parcels (`Port<SurfaceMsg>`);
//!   - the gas solver **consumes** them, advances a boundary-driven gas field in its
//!     own `CfdState`, samples the near-surface gas velocity at each surface parcel
//!     (`coupling::sample_gas_velocity`, verbatim) and forms the aerodynamic
//!     surface drag from the INDEPENDENT Schiller–Naumann `C_d(Re)` closure
//!     (`coupling::sphere_drag_force`, verbatim), then **exposes** that per-parcel
//!     drag (`Port<DragMsg>`);
//!   - the SPH bed **consumes** the drag and applies a first-principles surface
//!     entrainment law (drag vs submerged grain weight) to erode the free surface.
//!
//! This is the "couple any grass solver to any other" contract from
//! grass-solver-coupling-ergonomics (grass PR #11), exercised across two different
//! paradigms (Eulerian mesh gas ↔ Lagrangian granular continuum).
//!
//! ## Modelling context and controls (not external PSI validation)
//!
//! A surface grain is entrained when the aerodynamic drag overcomes the resisting
//! moment of its **submerged weight** — the Bagnold (1941) / Shields (1936)
//! incipient-motion criterion. Non-dimensionalised, the onset friction velocity is
//! `u*_t = A · sqrt((ρ_s−ρ_f)/ρ_f · g · d)` with the aerodynamic coefficient
//! `A ≈ 0.1` (Bagnold 1941; Iversen & White 1982). The following diagnostic checks
//! are assembled from Schiller–Naumann drag, bed friction, and a rough-wall log-law.
//! They do not supply a geometry-matched crater/erosion datum and therefore are not
//! acceptance evidence for the PSI goal:
//!   1. **Bagnold-A anchor.** The recovered `A_meas` (from the coupled seam's own
//!      onset slip, converted to `u*` by the textbook log-law) must land in the
//!      published context band around `A≈0.1`.
//!   2. **Reduced-gravity trend.** The onset slip must scale `u_gc ∝ g^{~1/2}` over
//!      Moon→Mars→Earth gravity (the PSI reduced-gravity behaviour). A **cohesive
//!      negative control** — resistance made *g-independent* — gives exponent ≈ 0
//!      and FAILS the band, proving the gate is not vacuous.
//!
//! ## The live coupled demonstration (Part C — the plume-surface interaction)
//!
//! With the baseline grain, the CFD inlet is stepped against the settled SPH bed
//! through the ports over a sweep of inlet strengths `U_peak = factor · u_gc`:
//!   - **Falsifiable threshold.** Above onset the surface **mobilises** (crater
//!     grows); below onset it stays **packed**.
//!   - **Location diagnostic.** The reported mean erosion offset is observational
//!     only. No Roberts location datum is digitized here, so it is not a gate.
//!   - **Negative control (falsifiable).** With the drag **port severed**, the SAME
//!     super-onset jet leaves the bed packed — the crater is caused by the coupling,
//!     not by numerics.
//!
//! ## Honest labelling
//!
//! The surface entrainment law (drag vs submerged weight) is the standard
//! Shields/Roberts closure, so "mobilises above / packed below onset" in Part C is a
//! **consistency** statement with Part B's onset. The reduced-gravity and
//! severed-port controls can falsify the implementation. None is a substitute for an external,
//! quantitative PSI comparison.
//!
//! ```text
//! cargo run --release -p cfd_ibm --example jet_crater_erosion -- \
//!     crates/cfd_ibm/examples/jet_crater_erosion.toml
//! ```
//!
//! References: R. A. Bagnold, *The Physics of Blown Sand and Desert Dunes*, Methuen
//! (1941); A. Shields, *Mitt. Preuss. Versuchsanst. Wasserbau Schiffbau* 26 (1936);
//! J. D. Iversen & B. R. White, *Sedimentology* 29:111 (1982); L. Roberts, "The
//! action of a hypersonic jet on a dust layer," IAS Paper 63-50 (1963);
//! L. Schiller & A. Naumann, *Z. Ver. Deut. Ing.* 77:318 (1935).

use std::any::TypeId;
use std::f64::consts::PI;

use cfd_boundary::{
    BoundaryPlugin, BoundaryRegistry, NoSlipWall, SubsonicInflow, SupersonicOutflow,
};
use cfd_eos::{Eos, EosResource, IdealGas, Viscosity};
use cfd_ibm::coupling::{self, cd_schiller_naumann, sphere_drag_force};
use cfd_solver::{
    CfdStatePlugin, FluxPlugin, IdealGasPlugin, IntegratorPlugin, SolverConfig, SolverPlugin,
};
use cfd_state::{CfdState, PrimVar};
use field_core::{
    FieldDefaultPlugins, FieldRegistry, MeshScheduleSet, UniformMesh, UniformMeshConfig, Vec3,
};
use grass_app::prelude::*;
use grass_multi::{consume_field, expose_field, tick_subapp, MultiAppExt, SubApps};
use serde::Deserialize;
use sph_core::prelude::*;

const R_GAS: f64 = 287.058; // matches IdealGas::air()

// ─── Declarative case ────────────────────────────────────────────────────────

#[derive(Deserialize, Default)]
struct GasCfg {
    rho: f64,
    p: f64,
    mu: f64,
}

#[derive(Deserialize, Default, Clone)]
struct BedCfg {
    rho_s: f64,
    rho_c: f64,
    mu_s: f64,
    mu_2: f64,
    i0: f64,
    bulk_modulus: f64,
    poisson: f64,
    grain_d: f64,
    restitution: f64,
    rest_density: f64,
    spacing: f64,
    x_lo: f64,
    x_hi: f64,
    y_lo: f64,
    y_hi: f64,
    z_lo: f64,
    z_hi: f64,
    floor_thickness: f64,
    settle_steps: usize,
    sph_dt: f64,
}

#[derive(Deserialize, Default)]
struct GridCfg {
    nx: usize,
    ny: usize,
    nz: usize,
    ng: usize,
    z_hi: f64,
}

#[derive(Deserialize, Default, Clone, Copy)]
struct JetCfg {
    x_center: f64,
    a: f64,
    surface_band: f64,
}

#[derive(Deserialize, Default)]
struct GravityCfg {
    gz: f64,
}

#[derive(Deserialize, Default)]
struct LogLawCfg {
    kappa: f64,
    roughness_ratio: f64,
}

#[derive(Deserialize, Default)]
struct ScalingCfg {
    g_list: Vec<f64>,
}

#[derive(Deserialize, Default)]
struct RunCfg {
    u_factors: Vec<f64>,
    dyn_steps: usize,
}

#[derive(Deserialize, Default)]
struct ValidationCfg {
    bagnold_a_ref: f64,
    bagnold_a_lo: f64,
    bagnold_a_hi: f64,
    bagnold_a_err_floor: f64,
    grav_exponent_lo: f64,
    grav_exponent_hi: f64,
    erode_min_hspeed: f64,
}

// ─── Independent surface-grain entrainment force balance (Bagnold/Shields) ─────
// The resisting force is the incipient-motion criterion: the aerodynamic drag on an
// exposed grain must overcome tan(φ)·(submerged weight), with tan(φ)=μ_s the bed's
// own repose friction (a material property, not a fitted constant). The DRAG is the
// INDEPENDENT Schiller–Naumann C_d(Re) sphere drag reused verbatim from cfd_ibm.

/// Onset near-surface gas slip `u_gc` at which an exposed surface grain is entrained.
/// `cohesive` replaces the submerged-weight resistance with a fixed (g-independent)
/// value — the reduced-gravity negative control.
fn entrainment_onset_u(
    d: f64,
    rho_s: f64,
    rho_f: f64,
    mu: f64,
    g: f64,
    mu_s: f64,
    cohesive_resist: Option<f64>,
) -> f64 {
    let r = 0.5 * d;
    let v = 4.0 / 3.0 * PI * r.powi(3);
    // Resisting force (horizontal drag needed to dislodge the grain).
    let resist = cohesive_resist.unwrap_or(mu_s * (rho_s - rho_f) * g * v);
    // Horizontal drag magnitude on a resting grain in slip `u` (monotone in u).
    let drag_mag = |u: f64| {
        let f = sphere_drag_force([u, 0.0, 0.0], [0.0; 3], r, rho_f, mu, cd_schiller_naumann);
        (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()
    };
    // Bisect drag_mag(u) = resist.
    let (mut lo, mut hi) = (1e-4, 500.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if drag_mag(mid) < resist {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Bagnold coefficient `A = u*/sqrt((ρ_s−ρ_f)/ρ_f · g · d)` recovered from the onset
/// slip `u_gc` via the textbook rough-wall log-law `u* = u_gc·κ/ln(z_ref/z0)`,
/// `z_ref=d`, `z0=d/roughness_ratio`. κ and z0/d are universal, not fitted.
fn bagnold_a(u_gc: f64, d: f64, rho_s: f64, rho_f: f64, g: f64, ll: &LogLawCfg) -> f64 {
    let u_star = u_gc * ll.kappa / ll.roughness_ratio.ln();
    u_star / (((rho_s - rho_f) / rho_f) * g * d).sqrt()
}

/// Least-squares slope of `ln(y)` vs `ln(x)` — the power-law exponent.
fn loglog_slope(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let lx: Vec<f64> = xs.iter().map(|v| v.ln()).collect();
    let ly: Vec<f64> = ys.iter().map(|v| v.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..xs.len() {
        num += (lx[i] - mx) * (ly[i] - my);
        den += (lx[i] - mx) * (lx[i] - mx);
    }
    num / den
}

// ─── The grass coupling contract: the two port message types ──────────────────
// These are the ONLY types shared between the gas and granular solvers.

/// SPH → gas: the free-surface parcels the gas may erode (index, position, velocity,
/// grain radius, and the entrainment-resisting force μ_s·(submerged weight)).
#[derive(Clone, Default)]
struct SurfaceMsg {
    idx: Vec<usize>,
    pos: Vec<Vec3>,
    vel: Vec<Vec3>,
    radius: Vec<f64>,
    resist: Vec<f64>,
}

/// gas → SPH: the aerodynamic drag on each exposed surface parcel, plus the sampled
/// near-surface gas speed there (diagnostic for the erosion-localisation measure).
#[derive(Clone, Default)]
struct DragMsg {
    idx: Vec<usize>,
    drag: Vec<Vec3>,
    u_local: Vec<f64>,
    x: Vec<f64>,
}

// ─── SPH sub-App (the granular free surface) ──────────────────────────────────

fn sph_config_toml(bed: &BedCfg, gz: f64, dt: f64) -> String {
    let floor_lo = bed.z_lo - bed.floor_thickness;
    let dom_z_hi = bed.z_hi + 6.0 * bed.spacing;
    let bin = 4.0 * bed.spacing;
    format!(
        "[comm]\nprocessors_x=1\nprocessors_y=1\nprocessors_z=1\n\
         [domain]\nx_low={xlo}\nx_high={xhi}\ny_low={ylo}\ny_high={yhi}\nz_low={fzlo}\nz_high={zhi}\n\
         boundary_x=\"fixed\"\nboundary_y=\"periodic\"\nboundary_z=\"fixed\"\n\
         [neighbor]\nnewton=false\nskin_fraction=1.2\nbin_size={bin}\n\
         [gravity]\ngx=0.0\ngy=0.0\ngz={gz}\n\
         [sph]\n[[sph.materials]]\nname=\"grain\"\nmu_s={mus}\nmu_2={mu2}\ni0={i0}\n\
         rho_s={rhos}\nrho_c={rhoc}\nbulk_modulus={bulk}\npoisson={pois}\nd={d}\nrestitution={rest}\n\
         [[sph.insert]]\nmaterial=\"grain\"\nregion_min=[{xlo},{ylo},{fzlo}]\nregion_max=[{xhi},{yhi},{zlo}]\n\
         spacing={sp}\nfrozen=true\n\
         [[sph.insert]]\nmaterial=\"grain\"\nregion_min=[{xlo},{ylo},{zlo}]\nregion_max=[{xhi},{yhi},{bzhi}]\n\
         spacing={sp}\nrest_density={rd}\n\
         [output]\ndir=\"/tmp/jet_crater_dump\"\n\
         [[run]]\nname=\"settle\"\nsteps=100000000\ndt={dt}\nthermo=0\n",
        xlo = bed.x_lo, xhi = bed.x_hi, ylo = bed.y_lo, yhi = bed.y_hi,
        zlo = bed.z_lo, zhi = dom_z_hi, fzlo = floor_lo, bzhi = bed.z_hi,
        bin = bin, gz = gz,
        mus = bed.mu_s, mu2 = bed.mu_2, i0 = bed.i0, rhos = bed.rho_s, rhoc = bed.rho_c,
        bulk = bed.bulk_modulus, pois = bed.poisson, d = bed.grain_d, rest = bed.restitution,
        sp = bed.spacing, rd = bed.rest_density, dt = dt,
    )
}

/// Static parameters the SPH coupling system needs (set on the SPH sub-App).
#[derive(Clone)]
struct SphCoupleParams {
    rho_s: f64,
    rho_f: f64,
    g: f64,
    mu_s: f64,
    surface_layer: f64, // free parcels within this of z_top are "surface"
    sever: bool,        // negative control: ignore the consumed drag (port severed)
}

/// Consumed drag from the gas port (index-aligned to this SPH app's atoms).
#[derive(Default)]
struct SphDrag {
    force: Vec<Vec3>, // per-atom horizontal drag, 0 where not eroded/eligible
}

/// The SPH free-surface parcels published to the gas port each step.
#[derive(Default, Clone)]
struct SphSurface(SurfaceMsg);

/// Live erosion diagnostic: the x-location and horizontal speed of the surface
/// parcels that were ACTUALLY entrained this step (drag beat the grain's resistance
/// AND the coupling was live). Empty ⇒ no erosion. Read back by `measure_crater`.
#[derive(Default)]
struct SphErosionDiag {
    x: Vec<f64>,
    hspeed: Vec<f64>,
}

/// SPH `Force` phase: (1) refresh the surface-parcel list for export, and (2) apply
/// the consumed gas drag to surface parcels **iff** it overcomes the grain's
/// entrainment resistance (Shields/Roberts surface law) — eroding the free surface.
fn sph_surface_and_erode(
    mut atoms: ResMut<Atom>,
    registry: Res<AtomDataRegistry>,
    params: Res<SphCoupleParams>,
    drag: Res<SphDrag>,
    mut surf: ResMut<SphSurface>,
    mut diag: ResMut<SphErosionDiag>,
) {
    let sph = registry.expect::<SphAtom>("sph_surface_and_erode");
    let n = atoms.nlocal as usize;
    // Current free-surface top.
    let mut z_top = f64::NEG_INFINITY;
    for i in 0..n {
        if sph.is_boundary[i] < 0.5 {
            z_top = z_top.max(atoms.pos[i][2]);
        }
    }
    let mut msg = SurfaceMsg::default();
    diag.x.clear();
    diag.hspeed.clear();
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        let is_surface = atoms.pos[i][2] >= z_top - params.surface_layer;
        if !is_surface {
            continue;
        }
        let v_solid = sph.particle_mass[i] / params.rho_s;
        let radius = (3.0 * v_solid / (4.0 * PI)).cbrt();
        let resist = params.mu_s * (params.rho_s - params.rho_f) * params.g * v_solid;
        // Apply the consumed gas drag with the entrainment gate (Shields/Roberts
        // surface law): a grain erodes only when the drag beats its resistance.
        if !params.sever {
            if let Some(f) = drag.force.get(i) {
                let fmag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
                if fmag > resist {
                    atoms.force[i][0] += f[0];
                    atoms.force[i][1] += f[1];
                    atoms.force[i][2] += f[2];
                    let v = atoms.vel[i];
                    diag.x.push(atoms.pos[i][0]);
                    diag.hspeed.push((v[0] * v[0] + v[1] * v[1]).sqrt());
                }
            }
        }
        msg.idx.push(i);
        msg.pos.push(atoms.pos[i]);
        msg.vel.push(atoms.vel[i]);
        msg.radius.push(radius);
        msg.resist.push(resist);
    }
    surf.0 = msg;
}

fn build_sph_app(bed: &BedCfg, gz: f64, dt: f64, params: SphCoupleParams) -> App {
    let toml = sph_config_toml(bed, gz, dt);
    let mut app = App::new();
    app.add_resource(grass_io::Config::from_str(&toml));
    app.add_resource(grass_io::Input {
        filename: String::from("jet_crater_sph"),
        output_dir: Some(String::from("/tmp/jet_crater_dump")),
    });
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(SphGravityPlugin);
    app.add_resource(params);
    app.add_resource(SphDrag::default());
    app.add_resource(SphSurface::default());
    app.add_resource(SphErosionDiag::default());
    app.add_update_system(sph_surface_and_erode, ParticleSimScheduleSet::Force);
    app
}

// ─── Gas sub-App (the boundary-driven impinging inflow) ────────────────────────

#[derive(Clone, Copy)]
struct GasProps {
    mu: f64,
    u_peak: f64,
}

/// The surface parcels the gas consumed from the SPH port this step.
#[derive(Default)]
struct GasConsumedSurface(SurfaceMsg);

/// The per-parcel drag the gas exposes back to the SPH port this step.
#[derive(Default, Clone)]
struct GasDrag(DragMsg);

/// Gas `Output`: sample the advanced CFD state and form per-parcel drag.  It never
/// overwrites `CfdState`; flux, boundary, CFL, and RK plugins own gas evolution.
fn gas_impose_and_drag(
    mesh: Res<UniformMesh>,
    reg: Res<FieldRegistry>,
    eos: Res<EosResource>,
    gas: Res<GasProps>,
    surf: Res<GasConsumedSurface>,
    mut out: ResMut<GasDrag>,
) {
    let eos: &dyn Eos = &*eos.0;
    let state = reg.expect_mut::<CfdState>("CfdState");
    let s = &surf.0;
    // Skip until the SPH has published a free surface.
    if s.pos.is_empty() {
        out.0 = DragMsg::default();
        return;
    }
    // Sample the evolved gas at each surface parcel and form the drag.
    let mut msg = DragMsg::default();
    for k in 0..s.idx.len() {
        let p = s.pos[k];
        let u_gas = coupling::sample_gas_velocity(&*mesh, &state, eos, p)
            .expect("surface parcel must lie inside the CFD mesh");
        let rho_f = coupling::sample_gas_density(&*mesh, &state, p)
            .expect("surface parcel must lie inside the CFD mesh");
        let f = sphere_drag_force(
            u_gas,
            s.vel[k],
            s.radius[k],
            rho_f,
            gas.mu,
            cd_schiller_naumann,
        );
        msg.idx.push(s.idx[k]);
        msg.drag.push(f);
        msg.u_local
            .push((u_gas[0] * u_gas[0] + u_gas[1] * u_gas[1] + u_gas[2] * u_gas[2]).sqrt());
        msg.x.push(p[0]);
    }
    out.0 = msg;
}

fn build_gas_app(gas: &GasCfg, mesh_cfg: UniformMeshConfig, props: GasProps) -> App {
    let (rho, p) = (gas.rho, gas.p);
    let t = p / (rho * R_GAS);
    let init = move |_x: Vec3| {
        let eos = IdealGas::air();
        // Start quiescent: the physical inlet below, not an interior write,
        // supplies the jet. The CFD update determines its route to the bed.
        eos.prim_to_cons(&PrimVar::new(rho, 0.0, 0.0, 0.0, p, t))
    };
    let bcs = BoundaryRegistry::default()
        .with_axis(0, SupersonicOutflow)
        .with_axis(1, SupersonicOutflow)
        .with(field_core::BoundarySide::ZLo, NoSlipWall)
        .with(
            field_core::BoundarySide::ZHi,
            SubsonicInflow {
                rho: gas.rho,
                u: 0.0,
                v: 0.0,
                w: -props.u_peak,
                t,
            },
        );
    let mut app = App::new();
    app.add_plugins(FieldDefaultPlugins { mesh: mesh_cfg })
        .add_plugins(CfdStatePlugin::new(init))
        .add_plugins(IdealGasPlugin)
        .add_plugins(BoundaryPlugin::<UniformMesh>::new(bcs))
        .add_plugins(FluxPlugin::<UniformMesh>::hllc())
        .add_plugins(IntegratorPlugin::rk3())
        .add_plugins(SolverPlugin::<UniformMesh>::new(SolverConfig {
            cfl: 0.25,
            muscl: false,
            viscous: true,
            ..SolverConfig::default()
        }));
    app.add_resource(Viscosity::Constant(gas.mu));
    app.add_resource(props);
    app.add_resource(GasConsumedSurface::default());
    app.add_resource(GasDrag::default());
    app.add_update_system(gas_impose_and_drag, MeshScheduleSet::Output);
    app
}

// ─── The coupled parent: two solvers wired ONLY through grass ports ───────────

#[derive(Debug, Clone, Copy)]
enum Phase {
    TickSph,     // SPH steps (+ refreshes its surface, applies last step's drag)
    ExposeSurf,  // SPH surface → Port<SurfaceMsg>
    ConsumeSurf, // Port<SurfaceMsg> → gas
    TickGas,     // gas imposes wall jet + forms drag
    ExposeDrag,  // gas drag → Port<DragMsg>
    ConsumeDrag, // Port<DragMsg> → SPH
}
impl ScheduleSet for Phase {
    fn to_index(&self) -> u32 {
        match self {
            Self::TickSph => 0,
            Self::ExposeSurf => 1,
            Self::ConsumeSurf => 2,
            Self::TickGas => 3,
            Self::ExposeDrag => 4,
            Self::ConsumeDrag => 5,
        }
    }
    fn name(&self) -> &'static str {
        match self {
            Self::TickSph => "TickSph",
            Self::ExposeSurf => "ExposeSurf",
            Self::ConsumeSurf => "ConsumeSurf",
            Self::TickGas => "TickGas",
            Self::ExposeDrag => "ExposeDrag",
            Self::ConsumeDrag => "ConsumeDrag",
        }
    }
}

struct CraterResult {
    /// Mean horizontal speed of the ENTRAINED (eroding) surface parcels.
    mean_eroding_hspeed: f64,
    /// Mean |x−x_center| of the entrained parcels (Roberts off-axis offset).
    erosion_offset: f64,
    n_surface: usize,
    n_eroding: usize,
}

/// One coupled run at wall-jet strength `u_peak`: build SPH+gas sub-Apps, settle the
/// bed, march the parent loop `steps` times through the ports, and measure the crater.
#[allow(clippy::too_many_arguments)]
fn run_coupled(
    gas: &GasCfg,
    bed: &BedCfg,
    grid: &GridCfg,
    jet: &JetCfg,
    g: f64,
    gz: f64,
    u_peak: f64,
    steps: usize,
    sever: bool,
) -> CraterResult {
    let dt = bed.sph_dt;
    let params = SphCoupleParams {
        rho_s: bed.rho_s,
        rho_f: gas.rho,
        g,
        mu_s: bed.mu_s,
        surface_layer: 1.2 * bed.spacing,
        sever,
    };
    let sph = build_sph_app(bed, gz, dt, params);
    let mesh_cfg = UniformMeshConfig {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        ng: grid.ng,
        bounds_lo: [bed.x_lo, bed.y_lo, bed.z_lo],
        bounds_hi: [bed.x_hi, bed.y_hi, grid.z_hi],
        y_edges: None,
        z_edges: None,
    };
    let props = GasProps { mu: gas.mu, u_peak };
    let gas_app = build_gas_app(gas, mesh_cfg, props);

    let mut parent = App::new();
    parent.add_subapp("sph", sph);
    parent.add_subapp("gas", gas_app);
    parent.add_port::<SurfaceMsg>();
    parent.add_port::<DragMsg>();

    parent.add_update_system(tick_subapp("sph", 1), Phase::TickSph);
    // SPH surface → gas, through Port<SurfaceMsg>.
    parent.add_update_system(
        expose_field::<SphSurface, SurfaceMsg>("sph", |s| s.0.clone()),
        Phase::ExposeSurf,
    );
    parent.add_update_system(
        consume_field::<GasConsumedSurface, SurfaceMsg>("gas", |c, m| c.0 = m.clone()),
        Phase::ConsumeSurf,
    );
    parent.add_update_system(tick_subapp("gas", 1), Phase::TickGas);
    // gas drag → SPH, through Port<DragMsg>. `apply` scatters the port drag into the
    // consumer's per-atom `SphDrag` (index-aligned); the SPH Force system gates it.
    parent.add_update_system(
        expose_field::<GasDrag, DragMsg>("gas", |g| g.0.clone()),
        Phase::ExposeDrag,
    );
    parent.add_update_system(
        consume_field::<SphDrag, DragMsg>("sph", scatter_drag),
        Phase::ConsumeDrag,
    );
    parent.prepare();

    // Settle the SPH bed before the jet is applied.
    {
        let cell = parent.get_mut_resource(TypeId::of::<SubApps>()).unwrap();
        let mut gd = cell.borrow_mut();
        let subs = gd.downcast_mut::<SubApps>().unwrap();
        for _ in 0..bed.settle_steps {
            subs.tick("sph");
        }
    }
    for _ in 0..steps {
        parent.run();
    }

    let result = measure_crater(&parent, jet);
    if let Some(cell) = parent.get_mut_resource(TypeId::of::<SubApps>()) {
        cell.borrow_mut()
            .downcast_mut::<SubApps>()
            .unwrap()
            .cleanup_all();
    }
    result
}

/// Consumer `apply`: scatter the port's per-parcel drag into a dense per-atom vector
/// sized to the SPH app so the `Force` system can index it by atom id.
fn scatter_drag(dst: &mut SphDrag, msg: &DragMsg) {
    let max_idx = msg.idx.iter().copied().max().unwrap_or(0);
    dst.force.clear();
    dst.force.resize(max_idx + 1, [0.0; 3]);
    for k in 0..msg.idx.len() {
        dst.force[msg.idx[k]] = msg.drag[k];
    }
}

/// Read the SPH sub-App's erosion diagnostic to score the crater.
fn measure_crater(parent: &App, jet: &JetCfg) -> CraterResult {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let sph = subs.find("sph").unwrap();
    let n_surface = {
        let surf_cell = sph
            .resource_cell(TypeId::of::<SphSurface>())
            .unwrap()
            .borrow();
        surf_cell.downcast_ref::<SphSurface>().unwrap().0.idx.len()
    };
    let diag_cell = sph
        .resource_cell(TypeId::of::<SphErosionDiag>())
        .unwrap()
        .borrow();
    let diag = diag_cell.downcast_ref::<SphErosionDiag>().unwrap();
    let n_eroding = diag.x.len();
    let mean_h = if n_eroding > 0 {
        diag.hspeed.iter().sum::<f64>() / n_eroding as f64
    } else {
        0.0
    };
    let offset = if n_eroding > 0 {
        diag.x.iter().map(|x| (x - jet.x_center).abs()).sum::<f64>() / n_eroding as f64
    } else {
        0.0
    };
    CraterResult {
        mean_eroding_hspeed: mean_h,
        erosion_offset: offset,
        n_surface,
        n_eroding,
    }
}

// ─── main ─────────────────────────────────────────────────────────────────────

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: jet_crater_erosion <case.toml>");
    let toml_src =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let cfg = grass_io::Config::from_str(&toml_src);

    let gas: GasCfg = cfg.section("gas");
    let bed: BedCfg = cfg.section("bed");
    let grid: GridCfg = cfg.section("grid");
    let jet: JetCfg = cfg.section("jet");
    let grav: GravityCfg = cfg.section("gravity");
    let ll: LogLawCfg = cfg.section("loglaw");
    let scaling: ScalingCfg = cfg.section("scaling");
    let run: RunCfg = cfg.section("run");
    let valid: ValidationCfg = cfg.section("validation");
    let g0 = grav.gz.abs();

    println!("# Jet–crater erosion: CFD gas wall-jet driving a granular-SPH free surface");
    println!("# COUPLING: grass exchange ports (grass_multi::Port / expose_field / consume_field)");
    println!("# CONTEXT: Bagnold (1941) threshold + Roberts (1963) wall-jet mechanism");
    println!(
        "# gas: rho_f={} mu={:.3e}   grain d={:.3e} m   rho_s={}   mu_s(tanφ)={}",
        gas.rho, gas.mu, bed.grain_d, bed.rho_s, bed.mu_s
    );

    // ── Part B1: independent Bagnold-A anchor at baseline ────────────────────────
    let u_gc = entrainment_onset_u(bed.grain_d, bed.rho_s, gas.rho, gas.mu, g0, bed.mu_s, None);
    let re_c = gas.rho * u_gc * bed.grain_d / gas.mu;
    let a_meas = bagnold_a(u_gc, bed.grain_d, bed.rho_s, gas.rho, g0, &ll);
    let a_err = (a_meas - valid.bagnold_a_ref).abs() / valid.bagnold_a_ref;
    println!("#");
    println!("# ── Part B: Bagnold/Shields mechanism diagnostics (not PSI validation) ──");
    println!(
        "# onset slip u_gc={u_gc:.4} m/s   Re_c={re_c:.0}   u*={:.4} m/s",
        u_gc * ll.kappa / ll.roughness_ratio.ln()
    );
    println!(
        "# Bagnold A_meas={a_meas:.4}  (band [{:.3},{:.3}], ref {:.2})  rel.err vs 0.10 = {:.1}% (floor {:.0}%)",
        valid.bagnold_a_lo, valid.bagnold_a_hi, valid.bagnold_a_ref, 100.0 * a_err, 100.0 * valid.bagnold_a_err_floor
    );

    // ── Part B2: reduced-gravity trend u_gc ∝ g^p + cohesive negative control ─────
    let mut gs = Vec::new();
    let mut us = Vec::new();
    let mut us_coh = Vec::new();
    // Cohesive control: freeze the resistance at its g0 value (g-independent).
    let coh = {
        let r = 0.5 * bed.grain_d;
        let v = 4.0 / 3.0 * PI * r.powi(3);
        bed.mu_s * (bed.rho_s - gas.rho) * g0 * v
    };
    for &gg in &scaling.g_list {
        gs.push(gg);
        us.push(entrainment_onset_u(
            bed.grain_d,
            bed.rho_s,
            gas.rho,
            gas.mu,
            gg,
            bed.mu_s,
            None,
        ));
        us_coh.push(entrainment_onset_u(
            bed.grain_d,
            bed.rho_s,
            gas.rho,
            gas.mu,
            gg,
            bed.mu_s,
            Some(coh),
        ));
    }
    let p_g = loglog_slope(&gs, &us);
    let p_g_coh = loglog_slope(&gs, &us_coh);
    println!("# reduced-gravity trend u_gc ∝ g^p:");
    for i in 0..gs.len() {
        println!(
            "#   g={:>5.2}  u_gc={:.4}   (cohesive control u_gc={:.4})",
            gs[i], us[i], us_coh[i]
        );
    }
    println!(
        "# exponent p={p_g:.3}  (band [{:.2},{:.2}], Bagnold ½)   cohesive-control p={p_g_coh:.3}  (must fail band)",
        valid.grav_exponent_lo, valid.grav_exponent_hi
    );

    // ── Part C: live coupled erosion through the grass ports ─────────────────────
    println!("#");
    println!("# ── Part C: live coupled wall-jet → SPH crater (through grass ports) ──");
    println!("#   U/u_gc   U_peak[m/s]   n_erode/n_surf   mean eroding |v_h|   offset/a   state");
    let mut below_ok = true;
    let mut below_seen = false;
    let mut above_ok = true;
    let mut above_seen = false;
    for &fac in &run.u_factors {
        let u_peak = fac * u_gc;
        let r = run_coupled(
            &gas,
            &bed,
            &grid,
            &jet,
            g0,
            grav.gz,
            u_peak,
            run.dyn_steps,
            false,
        );
        let eroding = r.n_eroding > 0 && r.mean_eroding_hspeed > valid.erode_min_hspeed;
        let off = if r.n_eroding > 0 {
            r.erosion_offset / jet.a
        } else {
            0.0
        };
        let state = if fac > 1.0 {
            above_seen = true;
            if !eroding {
                above_ok = false;
            }
            if eroding {
                "ERODES (exploratory)"
            } else {
                "no erosion (unexpected)"
            }
        } else {
            below_seen = true;
            if r.n_eroding > 0 {
                below_ok = false;
            }
            if r.n_eroding == 0 {
                "packed (bed intact)"
            } else {
                "eroded (unexpected)"
            }
        };
        println!(
            "  {fac:>6.2}   {u_peak:>10.4}   {:>6}/{:<6}   {:>18.4e}   {off:>8.3}   {state}",
            r.n_eroding, r.n_surface, r.mean_eroding_hspeed
        );
    }

    // Negative control: severed drag port at the strongest jet — must NOT erode.
    let strongest = run
        .u_factors
        .iter()
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let sev = run_coupled(
        &gas,
        &bed,
        &grid,
        &jet,
        g0,
        grav.gz,
        strongest * u_gc,
        run.dyn_steps,
        true,
    );
    let pass_sever = sev.n_eroding == 0;
    println!(
        "# severed-port control (U/u_gc={strongest:.2}): n_erode={} -> {}",
        sev.n_eroding,
        if pass_sever {
            "packed (no crater; the coupling is the cause)"
        } else {
            "ERODED WITHOUT GAS — control vacuous!"
        }
    );

    // ── Verdict ──────────────────────────────────────────────────────────────────
    let pass_a_band = a_meas >= valid.bagnold_a_lo && a_meas <= valid.bagnold_a_hi;
    let pass_a_indep = a_err > valid.bagnold_a_err_floor;
    let pass_grav = p_g >= valid.grav_exponent_lo && p_g <= valid.grav_exponent_hi;
    let pass_grav_neg = p_g_coh < valid.grav_exponent_lo; // cohesive control must fail the band
    let pass_below = below_seen && below_ok;
    let pass_above = above_seen && above_ok;

    println!("#");
    println!("# ── result ─────────────────────────────────────────────");
    println!(
        "# checks: A_band={pass_a_band} A_indep={pass_a_indep} grav_exp={pass_grav} grav_neg={pass_grav_neg} \
         packed_below={pass_below} erodes_above={pass_above} severed_neg={pass_sever}; erosion offset is diagnostic only"
    );
    if pass_a_band
        && pass_a_indep
        && pass_grav
        && pass_grav_neg
        && pass_below
        && pass_above
        && pass_sever
    {
        println!(
            "EXPLORATORY CONTROLS: PASS  (Bagnold/context and coupled-seam controls pass; \
             this is NOT external PSI validation: no digitized geometry-matched crater/erosion data are compared)",
        );
    } else {
        println!("EXPLORATORY CONTROLS: FAIL");
        std::process::exit(1);
    }
}
