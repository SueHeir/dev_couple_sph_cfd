//! **SPH--FIELD packed-bed force-transfer smoke case.**
//!
//! This program is an executable smoke case for a homogeneous packed bed. It imposes
//! an interstitial gas velocity; it does not solve a nozzle flow, advance an
//! impinging plume, form a crater, or predict erosion/ejecta. `U_mf` is reported
//! as a diagnostic for the force-transfer seam, not accepted as plume-surface
//! validation.
//!
//! ## What is new here (the goal): the granular phase is a μ(I) SPH CONTINUUM
//!
//! Unlike the merged `fluidized_bed_umf` capstone — where the bed is a packing of
//! discrete `soil::Atom` grains — here the granular phase is a **μ(I) elasto-
//! viscoplastic continuum** carried by the dev_soil_sph SPH tier (`sph_core`): a real
//! settling bed of continuum parcels, each carrying many grains. The parcels are
//! coupled to the gas (dev_field_efvm FIELD substrate, `CfdState`) through the **DEM–CFD
//! interphase seam reused VERBATIM** (`cfd_ibm::coupling`): the parcels' solid
//! volume is charged onto the mesh as a per-cell **solid volume fraction**, the seam
//! `drag_force_from_beta` / `apply_momentum_sink` kernels form the interphase force,
//! and the reaction feeds back to the parcels — i.e. **dev_sph supplies the per-cell
//! solid volume fraction + velocity in place of discrete particles**, which is
//! exactly the coupling this goal asks for. Each parcel carries a *solid* volume
//! `V_p = m/ρ_s` (charged onto the mesh so the deposit reproduces the bulk solid
//! fraction), while the GRAIN diameter `d` fed to the packed-bed closure is the
//! physical grain size from config — the two are distinct (a parcel is many grains).
//!
//! ## What this smoke case does and does not establish
//!
//! The live seam uses the MacDonald et al. (1979) packed-bed closure. It reports
//! Wen--Yu (1966) and an identical-seam discrete packing only as diagnostic
//! comparators. The latter shares the coupling implementation and is therefore a
//! consistency check, not independent evidence. The executable fault controls
//! demonstrate sensitivity to two implementation mistakes; they are not an
//! experimental validation. This program deliberately has no numerical pass band
//! or scientific verdict: its successful exit means only that the configured
//! cross-substrate execution completed.
//!
//! A plume/crater claim is blocked pending the external, held-out, matched-observable
//! protocol in `EXTERNAL_VALIDATION.md`. That protocol requires an adversarial
//! wrong-coupling comparison and cannot be satisfied by this program's output.
//!
//! ## The total fluid force — drag AND the ∇P (generalized-buoyancy) force
//!
//! A parcel in a bed feels two fluid forces, and getting `U_mf` right needs **both**
//! (this is the physics the PR #5 first draft got wrong — it applied drag only, so
//! its coupled bed never lifted). (1) interphase **drag** `F = β V_p /(1−ε)·u_rel`
//! (Σ = ε·(dP/L)·V_bed), and (2) the **pressure-gradient / generalized-buoyancy
//! force** `+V_p β u_rel/ε` (Σ = (1−ε)·(dP/L)·V_bed). Only their sum `= (dP/L)·V_bed`
//! balances the buoyant bed weight at the correct velocity.
//!
//! ## Imposed flow
//!
//! As in `fluidized_bed_umf`, the interstitial gas velocity `u_g = U/ε` is **imposed**
//! in the coupled cells rather than obtained from a compressible flow solve: a
//! porosity-weighted momentum solve through a dense bed is the resolved-track story,
//! out of scope for the unresolved seam, and the near-incompressible transient is far
//! slower than the bed response. Everything case-specific is TOML.
//!
//! ```text
//! cargo run --release --example packed_bed_seam -- examples/packed_bed_seam/config.toml
//! ```
//!
//! References: S. Ergun, *Chem. Eng. Prog.* 48(2):89 (1952); C. Y. Wen & Y. H. Yu,
//! *AIChE J.* 12:610 (1966); I. F. MacDonald et al., *Ind. Eng. Chem. Fundam.*
//! 18(3):199 (1979); D. Gidaspow, *Multiphase Flow and Fluidization* (1994).

use std::f64::consts::PI;

use field_core::{UniformMesh, UniformMeshConfig, Vec3};
use grass_multi::MultiAppExt;
use serde::Deserialize;
use sph_cfd::prelude::*;
use sph_core::prelude::*;

// ─── Declarative case ────────────────────────────────────────────────────────

#[derive(Deserialize, Default)]
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
struct DemCfg {
    solid_fraction: f64,
    ncx: usize,
    ncy: usize,
    ncz: usize,
}

#[derive(Deserialize, Default)]
struct RunCfg {
    /// Superficial-velocity factors (× U_mf) for the live dynamical fluidization
    /// sweep — must bracket 1.0 (some below onset, some above).
    dyn_factors: Vec<f64>,
    /// Coupled steps integrated at each dynamical sweep point.
    dyn_steps: usize,
}

/// This example has a deliberately narrow forcing model.  Keeping that choice
/// in the case file, and rejecting anything else, prevents a renamed or edited
/// case from being mistaken for an advancing nozzle/plume calculation.
#[derive(Deserialize, Default)]
struct FlowCfg {
    model: String,
}

fn require_imposed_uniform_flow(flow: &FlowCfg) {
    const SUPPORTED: &str = "imposed_uniform_interstitial_velocity";
    assert_eq!(
        flow.model, SUPPORTED,
        "packed_bed_seam supports only \
         `{SUPPORTED}`, not a nozzle, plume, crater, or erosion model"
    );
}

// ─── Part A: the settled μ(I)-continuum SPH bed ───────────────────────────────

fn sph_config_toml(bed: &BedCfg, grav: &GravityCfg, dt: f64) -> String {
    let floor_lo = bed.z_lo - bed.floor_thickness;
    let dom_z_hi = bed.z_hi + 4.0 * bed.spacing;
    let bin = 4.0 * bed.spacing;
    format!(
        "[comm]\nprocessors_x=1\nprocessors_y=1\nprocessors_z=1\n\
         [domain]\nx_low={xlo}\nx_high={xhi}\ny_low={ylo}\ny_high={yhi}\nz_low={fzlo}\nz_high={zhi}\n\
         boundary_x=\"periodic\"\nboundary_y=\"periodic\"\nboundary_z=\"fixed\"\n\
         [neighbor]\nnewton=false\nskin_fraction=1.2\nbin_size={bin}\n\
         [gravity]\ngx=0.0\ngy=0.0\ngz={gz}\n\
         [sph]\n[[sph.materials]]\nname=\"grain\"\nmu_s={mus}\nmu_2={mu2}\ni0={i0}\n\
         rho_s={rhos}\nrho_c={rhoc}\nbulk_modulus={bulk}\npoisson={pois}\nd={d}\nrestitution={rest}\n\
         [[sph.insert]]\nmaterial=\"grain\"\nregion_min=[{xlo},{ylo},{fzlo}]\nregion_max=[{xhi},{yhi},{zlo}]\n\
         spacing={sp}\nfrozen=true\n\
         [[sph.insert]]\nmaterial=\"grain\"\nregion_min=[{xlo},{ylo},{zlo}]\nregion_max=[{xhi},{yhi},{bzhi}]\n\
         spacing={sp}\nrest_density={rd}\n\
         [output]\ndir=\"/tmp/packed_bed_seam_dump\"\n\
         [[run]]\nname=\"settle\"\nsteps=100000000\ndt={dt}\n",
        xlo = bed.x_lo, xhi = bed.x_hi, ylo = bed.y_lo, yhi = bed.y_hi,
        zlo = bed.z_lo, zhi = dom_z_hi, fzlo = floor_lo, bzhi = bed.z_hi,
        bin = bin, gz = grav.gz,
        mus = bed.mu_s, mu2 = bed.mu_2, i0 = bed.i0, rhos = bed.rho_s, rhoc = bed.rho_c,
        bulk = bed.bulk_modulus, pois = bed.poisson, d = bed.grain_d, rest = bed.restitution,
        sp = bed.spacing, rd = bed.rest_density, dt = dt,
    )
}

fn build_sph_app(bed: &BedCfg, grav: &GravityCfg, dt: f64) -> App {
    let toml = sph_config_toml(bed, grav, dt);
    let mut app = App::new();
    app.add_resource(grass_io::Config::from_str(&toml));
    app.add_resource(grass_io::Input {
        filename: String::from("packed_bed_seam_sph"),
        output_dir: Some(String::from("/tmp/packed_bed_seam_dump")),
    });
    app.add_plugins(CorePlugins)
        .add_plugins(SphDefaultPlugins)
        .add_plugins(SphGravityPlugin);
    app.add_resource(FluidForces::default());
    app.add_update_system(sph_fluid_force, ParticleSimScheduleSet::Force);
    app
}

// ─── Discrete DEM FCC reference packing (the resolved DEM–CFD cross-reference) ─

fn fcc_packing(nc: [usize; 3], a: f64) -> (Vec<Vec3>, Vec3) {
    let s = 0.25;
    let basis = [
        [s, s, s],
        [s, 0.5 + s, 0.5 + s],
        [0.5 + s, s, 0.5 + s],
        [0.5 + s, 0.5 + s, s],
    ];
    let mut pos = Vec::with_capacity(4 * nc[0] * nc[1] * nc[2]);
    for i in 0..nc[0] {
        for j in 0..nc[1] {
            for k in 0..nc[2] {
                for b in &basis {
                    pos.push([
                        (i as f64 + b[0]) * a,
                        (j as f64 + b[1]) * a,
                        (k as f64 + b[2]) * a,
                    ]);
                }
            }
        }
    }
    (pos, [nc[0] as f64 * a, nc[1] as f64 * a, nc[2] as f64 * a])
}

/// Build the discrete DEM reference bed as seam `Parcel`s (one grain = one parcel:
/// radius = grain radius, v_solid = grain volume) and its geometric voidage.
fn dem_parcels(dem: &DemCfg, d_grain: f64) -> (Vec<Parcel>, f64, Vec3) {
    let radius = 0.5 * d_grain;
    let v_solid = 4.0 / 3.0 * PI * radius.powi(3);
    let a = d_grain * (2.0 * PI / (3.0 * dem.solid_fraction)).cbrt();
    let (pos, bounds) = fcc_packing([dem.ncx, dem.ncy, dem.ncz], a);
    let n = pos.len();
    let v_bed = bounds[0] * bounds[1] * bounds[2];
    let eps = 1.0 - n as f64 * v_solid / v_bed;
    let parcels = pos
        .into_iter()
        .map(|center| Parcel { center, v_solid })
        .collect();
    (parcels, eps, bounds)
}

// ─── Gas mesh for the live coupled run ────────────────────────────────────────

fn build_gas_mesh(bed: &BedCfg, grid: &GridCfg) -> UniformMeshConfig {
    UniformMeshConfig {
        nx: grid.nx,
        ny: grid.ny,
        nz: grid.nz,
        ng: grid.ng,
        bounds_lo: [bed.x_lo, bed.y_lo, bed.z_lo],
        bounds_hi: [bed.x_hi, bed.y_hi, grid.z_hi],
        y_edges: None,
        z_edges: None,
    }
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: packed_bed_seam <case.toml>");
    let toml_src =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let cfg = grass_io::Config::from_str(&toml_src);

    let gas: GasCfg = cfg.section("gas");
    let bed: BedCfg = cfg.section("bed");
    let grid: GridCfg = cfg.section("grid");
    let grav: GravityCfg = cfg.section("gravity");
    let dem: DemCfg = cfg.section("dem");
    let run: RunCfg = cfg.section("run");
    let flow: FlowCfg = cfg.section("flow");
    require_imposed_uniform_flow(&flow);
    let g = grav.gz.abs();

    // ── Part A: settle the μ(I) continuum bed, measure ε_bed + settled contact p ──
    let mut sph = build_sph_app(&bed, &grav, bed.sph_dt);
    sph.prepare();
    for _ in 0..bed.settle_steps {
        sph.run();
    }
    let (parcels, eps_bed, z_top, z_bot, p_settled) = read_sph_bed(&sph, bed.rho_s);
    let n_parcels = parcels.len();

    // Deposit the settled continuum's solid volume onto a bed-scale mesh → per-cell
    // ε field (dev_sph supplying per-cell solid volume fraction to the gas seam). The
    // bulk column height is n_layers·spacing (each of the n_layers parcel layers
    // occupies one `spacing`), centred on the parcels so ε_cell ≈ ε_bed.
    let n_layers = ((z_top - z_bot) / bed.spacing).round() + 1.0;
    let z_height = n_layers * bed.spacing;
    let z_bottom = z_bot - 0.5 * bed.spacing;
    let dep_cfg = build_deposit_mesh(
        bed.x_lo,
        bed.x_hi,
        bed.y_lo,
        bed.y_hi,
        z_bottom,
        z_height,
        2.0 * bed.spacing,
    );
    let dep_mesh = UniformMesh::from_config(&dep_cfg);
    let (eps_field, cell_of_parcel) = deposit_bed_void_fraction(&dep_mesh, &parcels);
    let dep_err_sph = deposit_worst_err(&dep_mesh, &eps_field, eps_bed);

    // ── Part B: cross-method U_mf — SPH continuum vs discrete DEM, same live seam ──
    let mode = SeamMode::default();
    let umf_sph = measure_umf(
        &parcels,
        &eps_field,
        &cell_of_parcel,
        eps_bed,
        gas.rho,
        bed.rho_s,
        gas.mu,
        bed.grain_d,
        g,
        mode,
    );

    // Discrete DEM reference bed (the resolved DEM–CFD fluidization reference), same seam.
    let (dem_pv, eps_dem, dem_bounds) = dem_parcels(&dem, bed.grain_d);
    let dem_cell = (dem_bounds[0] / dem.ncx as f64).max(bed.grain_d);
    // FCC fills [0, bounds] exactly, so the bulk column height is bounds_z (z_bottom=0).
    let dem_dep_cfg = build_deposit_mesh(
        0.0,
        dem_bounds[0],
        0.0,
        dem_bounds[1],
        0.0,
        dem_bounds[2],
        dem_cell,
    );
    let dem_mesh = UniformMesh::from_config(&dem_dep_cfg);
    let (dem_eps_field, dem_cop) = deposit_bed_void_fraction(&dem_mesh, &dem_pv);
    let dep_err_dem = deposit_worst_err(&dem_mesh, &dem_eps_field, eps_dem);
    let umf_dem = measure_umf(
        &dem_pv,
        &dem_eps_field,
        &dem_cop,
        eps_dem,
        gas.rho,
        bed.rho_s,
        gas.mu,
        bed.grain_d,
        g,
        mode,
    );

    // Independent reference + analytic brackets.
    let (umf_wy, ar, re_mf) = u_mf_wen_yu(gas.rho, bed.rho_s, g, bed.grain_d, gas.mu);
    let umf_erg = u_mf_balance(
        150.0,
        1.75,
        eps_bed,
        gas.rho,
        bed.rho_s,
        g,
        bed.grain_d,
        gas.mu,
    );
    let umf_mac = u_mf_balance(
        180.0,
        1.8,
        eps_bed,
        gas.rho,
        bed.rho_s,
        g,
        bed.grain_d,
        gas.mu,
    );

    let err_dem = (umf_sph - umf_dem).abs() / umf_dem;
    let err_wy = (umf_sph - umf_wy).abs() / umf_wy;

    // Executable sensitivity probes. They are reported, never used as an
    // acceptance criterion: these closures and this forcing are not a held-out
    // plume-surface experiment.
    let umf_nopg = measure_umf(
        &parcels,
        &eps_field,
        &cell_of_parcel,
        eps_bed,
        gas.rho,
        bed.rho_s,
        gas.mu,
        bed.grain_d,
        g,
        SeamMode {
            omit_pressure_grad: true,
            ..mode
        },
    );
    let umf_epsbug = measure_umf(
        &parcels,
        &eps_field,
        &cell_of_parcel,
        eps_bed,
        gas.rho,
        bed.rho_s,
        gas.mu,
        bed.grain_d,
        g,
        SeamMode {
            corrupt_eps_power: true,
            ..mode
        },
    );
    println!("# SPH--FIELD packed-bed force-transfer smoke case (U_mf diagnostic)");
    println!("# MEASURED force: INDEPENDENT MacDonald et al. (1979) 180/1.8 closure via the seam (drag + grad-P)");
    println!(
        "# COMPARATORS:    discrete FCC uses the same seam; Wen & Yu is a packed-bed correlation"
    );
    println!(
        "# gas: rho={} mu={:.3e} p={:.3e}   grain d={:.3e} m   rho_s={} rho_f={}",
        gas.rho, gas.mu, gas.p, bed.grain_d, bed.rho_s, gas.rho
    );
    println!("# SPH bed: {n_parcels} free parcels   spacing={:.3e} m   z_top={:.4} m   settled contact p={:.3e} Pa", bed.spacing, z_top, p_settled);
    println!("# eps_bed(SPH continuum, settled)={eps_bed:.4}   eps_dem(FCC)={eps_dem:.4}   Ar={ar:.1}   Re_mf(WenYu)={re_mf:.3}");
    println!("#");
    println!("# ── minimum fluidization velocity U_mf [m/s] ───────────────────");
    println!("# SPH continuum (seam, MacDonald, drag+gradP) : {umf_sph:.5}   <- the coupled-continuum measurement");
    println!("# Wen & Yu 1966 correlation (packed-bed comparator) : {umf_wy:.5}   diagnostic difference {:.2}%", 100.0 * err_wy);
    println!("# DEM discrete  (same seam, MacDonald, drag+gradP) : {umf_dem:.5}   same-seam difference {:.2}%", 100.0 * err_dem);
    println!("# analytic brackets: Ergun(150/1.75)={umf_erg:.5} ({:+.2}% vs WenYu)  MacDonald(180/1.8)={umf_mac:.5} ({:+.2}% vs WenYu)", 100.0 * (umf_erg / umf_wy - 1.0), 100.0 * (umf_mac / umf_wy - 1.0));
    println!("# diagnostic spread vs WenYu {:.2}% (MacDonald and Wen--Yu are distinct packed-bed closures)", 100.0 * err_wy);
    println!("# sensitivity probes: omit-gradP {umf_nopg:.4} ({:+.1}%)  eps-power-bug {umf_epsbug:.4} ({:+.1}%)",
        100.0 * (umf_nopg / umf_wy - 1.0), 100.0 * (umf_epsbug / umf_wy - 1.0));
    println!(
        "# deposited-voidage diagnostic: worst |eps_cell-eps_bed|/eps_bed  SPH {:.2}%  DEM {:.2}%",
        100.0 * dep_err_sph,
        100.0 * dep_err_dem
    );

    // ── Part C: live coupled dynamical fluidization sweep (SPH continuum + gas) ──
    println!("#");
    println!("# ── live coupled dynamical fluidization (SPH continuum stepped in the seam) ──");
    println!("#   U/U_mf     U [m/s]     mean v_z [m/s]    contact p [Pa]   p/p_settled");
    let mut worst_mom = 0.0f64;
    for &fac in &run.dyn_factors {
        let u = fac * umf_sph;
        let dynr = run_coupled(&gas, &bed, &grid, &grav, eps_bed, u, run.dyn_steps);
        worst_mom = worst_mom.max(dynr.mom_err);
        println!(
            "  {fac:>7.2}   {u:>9.4}   {:>+13.4e}   {:>13.4e}   {:>10.4}",
            dynr.mean_vz,
            dynr.mean_p,
            dynr.mean_p / p_settled.max(1e-30)
        );
    }

    println!("# two-way momentum-exchange residual (diagnostic): {worst_mom:.2e}");
    println!("SMOKE COMPLETED: configured cross-substrate run executed; not a physical validation or plume-surface prediction.");
}

/// One coupled run at superficial `u_super`: build SPH+CFD sub-Apps, settle-prime the
/// SPH bed, march the parent loop `steps` times, and report mean free-parcel v_z,
/// mean contact pressure, and the worst two-way momentum-conservation error.
#[allow(clippy::too_many_arguments)]
fn run_coupled(
    gas: &GasCfg,
    bed: &BedCfg,
    grid: &GridCfg,
    grav: &GravityCfg,
    eps_bed: f64,
    u_super: f64,
    steps: usize,
) -> CoupledRunResult {
    let dt = bed.sph_dt;
    let props = CoupledGasProps {
        mu: gas.mu,
        rho: gas.rho,
        d_grain: bed.grain_d,
        eps_bed,
        dt,
        g: grav.gz.abs(),
        u_super,
    };
    let sph = build_sph_app(bed, grav, dt);
    let mesh_cfg = build_gas_mesh(bed, grid);
    let cfd = build_coupled_cfd(gas, mesh_cfg, props);

    let mut parent = App::new();
    parent.add_subapp("sph", sph);
    parent.add_subapp("cfd", cfd);
    add_standard_coupled_schedule(&mut parent);
    parent.prepare();

    prime_sph_and_spec(&mut parent, bed.settle_steps, bed.rho_s);
    for _ in 0..steps {
        parent.run();
    }
    let (mean_vz, mean_p) = sph_bed_state(&parent);
    let mom_err = coupled_diag(&parent).mom_err;
    cleanup_subapps(&mut parent);
    CoupledRunResult {
        mean_vz,
        mean_p,
        mom_err,
    }
}
