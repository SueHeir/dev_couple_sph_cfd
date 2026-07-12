use std::any::TypeId;
use std::f64::consts::PI;

use cfd_eos::{Eos, EosResource, IdealGas, Viscosity};
use cfd_ibm::coupling::{
    self, drag_force_from_beta, InterphaseForces, ParticleKinematics, ParticleSet,
};
use cfd_solver::{CfdStatePlugin, IdealGasPlugin};
use cfd_state::{CfdState, PrimVar};
use field_core::{
    FieldDefaultPlugins, FieldRegistry, FvMesh, MeshScheduleSet, UniformMesh, UniformMeshConfig,
    Vec3,
};
use grass_app::prelude::*;
use grass_multi::{tick_subapp, Multi, SubApps};
use sph_core::prelude::*;

use crate::bed::Parcel;
use crate::config::GasCfg;
use crate::drag::macdonald_beta;

pub const R_GAS: f64 = 287.058;

/// Per-parcel fluid load handed into the SPH sub-App, in `Atom` index order.
#[derive(Default)]
pub struct FluidForces {
    pub f: Vec<Vec3>,
}

/// `Force` phase on the SPH sub-App: add the seam's fluid load to each free parcel.
pub fn sph_fluid_force(
    mut atoms: ResMut<Atom>,
    ff: Res<FluidForces>,
    registry: Res<AtomDataRegistry>,
) {
    let sph = registry.expect::<SphAtom>("sph_fluid_force");
    let n = atoms.nlocal as usize;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        if let Some(f) = ff.f.get(i) {
            atoms.force[i][0] += f[0];
            atoms.force[i][1] += f[1];
            atoms.force[i][2] += f[2];
        }
    }
}

/// Read free SPH parcels as coupling parcels, plus mean voidage, bed top/bottom,
/// and mean grain-contact pressure.
pub fn read_sph_bed(app: &App, rho_s: f64) -> (Vec<Parcel>, f64, f64, f64, f64) {
    let atoms = app.get_resource_ref::<Atom>().expect("Atom");
    let registry = app
        .get_resource_ref::<AtomDataRegistry>()
        .expect("registry");
    let sph = registry.expect::<SphAtom>("read_sph_bed");
    let n = atoms.nlocal as usize;
    let mut parcels = Vec::new();
    let (mut sum_phi, mut sum_p) = (0.0, 0.0);
    let mut nfree = 0usize;
    let mut z_top = f64::NEG_INFINITY;
    let mut z_bot = f64::INFINITY;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        let v_solid = sph.particle_mass[i] / rho_s;
        parcels.push(Parcel {
            center: [atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2]],
            v_solid,
        });
        sum_phi += sph.density[i] / rho_s;
        sum_p += sph.pressure[i];
        nfree += 1;
        z_top = z_top.max(atoms.pos[i][2]);
        z_bot = z_bot.min(atoms.pos[i][2]);
    }
    let nf = nfree.max(1) as f64;
    (parcels, 1.0 - sum_phi / nf, z_top, z_bot, sum_p / nf)
}

#[derive(Default)]
pub struct ParcelSpec {
    pub radius: Vec<f64>,
    pub v_solid: Vec<f64>,
}

#[derive(Clone, Copy)]
pub struct CoupledGasProps {
    pub mu: f64,
    pub rho: f64,
    pub d_grain: f64,
    pub eps_bed: f64,
    pub dt: f64,
    pub g: f64,
    pub u_super: f64,
}

#[derive(Clone, Copy, Default)]
pub struct CoupledSeamDiag {
    pub mom_err: f64,
}

#[derive(Debug, Clone, Copy)]
pub enum Phase {
    Export,
    TickCfd,
    Import,
    TickSph,
}

impl ScheduleSet for Phase {
    fn to_index(&self) -> u32 {
        match self {
            Self::Export => 0,
            Self::TickCfd => 1,
            Self::Import => 2,
            Self::TickSph => 3,
        }
    }

    fn name(&self) -> &'static str {
        match self {
            Self::Export => "Export",
            Self::TickCfd => "TickCfd",
            Self::Import => "Import",
            Self::TickSph => "TickSph",
        }
    }
}

/// CFD sub-App `Output` phase: impose interstitial gas velocity, evaluate the
/// MacDonald seam force (drag + pressure-gradient + buoyancy), and apply the
/// equal-and-opposite drag momentum sink.
pub fn coupled_seam_system(
    mesh: Res<UniformMesh>,
    reg: Res<FieldRegistry>,
    eos: Res<EosResource>,
    gas: Res<CoupledGasProps>,
    spec: Res<ParcelSpec>,
    pset: Res<ParticleSet>,
    mut forces: ResMut<InterphaseForces>,
    mut diag: ResMut<CoupledSeamDiag>,
) {
    let eos: &dyn Eos = &*eos.0;
    let mut state = reg.expect_mut::<CfdState>("CfdState");
    let parts = &pset.particles;
    forces.reset(parts.len());
    if parts.is_empty() {
        return;
    }
    let eps = gas.eps_bed;
    let u_g_imp = gas.u_super / eps;
    for c in 0..mesh.n_cells_total() {
        if !mesh.is_local_cell(c) {
            continue;
        }
        let rho = state.u[c].rho;
        state.u[c].rho_u = 0.0;
        state.u[c].rho_v = 0.0;
        state.u[c].rho_w = rho * u_g_imp;
    }
    let mut drag_on_particle = vec![[0.0f64; 3]; parts.len()];
    for (i, p) in parts.iter().enumerate() {
        let u_gas = coupling::sample_gas_velocity(&*mesh, &state, eos, p.center)
            .unwrap_or([0.0, 0.0, u_g_imp]);
        let rho_f = coupling::sample_gas_density(&*mesh, &state, p.center).unwrap_or(gas.rho);
        let rel = [
            u_gas[0] - p.velocity[0],
            u_gas[1] - p.velocity[1],
            u_gas[2] - p.velocity[2],
        ];
        let rel_speed = (rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]).sqrt();
        let beta = macdonald_beta(eps, rho_f, gas.mu, gas.d_grain, rel_speed);
        let v_solid = spec.v_solid.get(i).copied().unwrap_or(0.0);
        let drag = drag_force_from_beta(beta, v_solid, eps, rel);
        let pg = v_solid * beta / eps;
        let buoy_z = rho_f * v_solid * gas.g;
        forces.force[i] = [
            drag[0] + pg * rel[0],
            drag[1] + pg * rel[1],
            drag[2] + pg * rel[2] + buoy_z,
        ];
        drag_on_particle[i] = drag;
    }

    let mut m0 = [0.0f64; 3];
    for c in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(c) {
            let v = mesh.cell_volume(c);
            m0[0] += state.u[c].rho_u * v;
            m0[1] += state.u[c].rho_v * v;
            m0[2] += state.u[c].rho_w * v;
        }
    }
    coupling::apply_momentum_sink(&*mesh, &mut state, parts, &drag_on_particle, gas.dt);
    let mut m1 = [0.0f64; 3];
    for c in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(c) {
            let v = mesh.cell_volume(c);
            m1[0] += state.u[c].rho_u * v;
            m1[1] += state.u[c].rho_v * v;
            m1[2] += state.u[c].rho_w * v;
        }
    }
    let (mut dn, mut sc) = (0.0f64, 0.0f64);
    for k in 0..3 {
        let dm = m1[k] - m0[k];
        let imp = -drag_on_particle.iter().map(|f| f[k]).sum::<f64>() * gas.dt;
        dn += (dm - imp) * (dm - imp);
        sc += imp * imp;
    }
    diag.mom_err = dn.sqrt() / sc.sqrt().max(1e-30);
}

pub fn build_coupled_cfd(gas: &GasCfg, mesh_cfg: UniformMeshConfig, props: CoupledGasProps) -> App {
    let (rho, p) = (gas.rho, gas.p);
    let t = p / (rho * R_GAS);
    let init = move |_x: Vec3| {
        let eos = IdealGas::air();
        eos.prim_to_cons(&PrimVar::new(rho, 0.0, 0.0, 0.0, p, t))
    };
    let mut app = App::new();
    app.add_plugins(FieldDefaultPlugins { mesh: mesh_cfg })
        .add_plugins(CfdStatePlugin::new(init))
        .add_plugins(IdealGasPlugin);
    app.add_resource(Viscosity::Constant(gas.mu));
    app.add_resource(props);
    app.add_resource(ParcelSpec::default());
    app.add_resource(ParticleSet::default());
    app.add_resource(InterphaseForces::default());
    app.add_resource(CoupledSeamDiag::default());
    app.add_update_system(coupled_seam_system, MeshScheduleSet::Output);
    app
}

pub fn export_kinematics(world: Multi) {
    let atoms = world.expect_read::<Atom>("sph");
    let n = atoms.nlocal as usize;
    let radii = {
        let spec = world.expect_read::<ParcelSpec>("cfd");
        spec.radius.clone()
    };
    let mut set = world.expect_write::<ParticleSet>("cfd");
    set.particles.clear();
    for i in 0..n {
        set.particles.push(ParticleKinematics {
            center: [atoms.pos[i][0], atoms.pos[i][1], atoms.pos[i][2]],
            velocity: [atoms.vel[i][0], atoms.vel[i][1], atoms.vel[i][2]],
            radius: radii.get(i).copied().unwrap_or(0.0),
        });
    }
}

pub fn import_forces(world: Multi) {
    let f = {
        let forces = world.expect_read::<InterphaseForces>("cfd");
        forces.force.clone()
    };
    world.expect_write::<FluidForces>("sph").f = f;
}

pub fn add_standard_coupled_schedule(parent: &mut App) {
    parent.add_update_system(export_kinematics, Phase::Export);
    parent.add_update_system(tick_subapp("cfd", 1), Phase::TickCfd);
    parent.add_update_system(import_forces, Phase::Import);
    parent.add_update_system(tick_subapp("sph", 1), Phase::TickSph);
}

/// Prime the SPH sub-App, then store constant per-parcel radii and solid volumes
/// as `ParcelSpec` on the CFD sub-App.
pub fn prime_sph_and_spec(parent: &mut App, settle_steps: usize, rho_s: f64) {
    {
        let cell = parent.get_mut_resource(TypeId::of::<SubApps>()).unwrap();
        let mut gd = cell.borrow_mut();
        let subs = gd.downcast_mut::<SubApps>().unwrap();
        for _ in 0..settle_steps {
            subs.tick("sph");
        }
    }
    let (radius, v_solid) = {
        let subs = parent.get_resource_ref::<SubApps>().unwrap();
        let sph = subs.find("sph").unwrap();
        let atom_cell = sph.resource_cell(TypeId::of::<Atom>()).unwrap().borrow();
        let reg_cell = sph
            .resource_cell(TypeId::of::<AtomDataRegistry>())
            .unwrap()
            .borrow();
        let atoms = atom_cell.downcast_ref::<Atom>().unwrap();
        let registry = reg_cell.downcast_ref::<AtomDataRegistry>().unwrap();
        let sph = registry.expect::<SphAtom>("prime_sph_and_spec");
        let n = atoms.nlocal as usize;
        let (mut radius, mut v_solid) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for i in 0..n {
            let vs = sph.particle_mass[i] / rho_s;
            if sph.is_boundary[i] > 0.5 {
                radius.push(0.0);
                v_solid.push(0.0);
            } else {
                radius.push((3.0 * vs / (4.0 * PI)).cbrt());
                v_solid.push(vs);
            }
        }
        (radius, v_solid)
    };
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let cfd = subs.find("cfd").unwrap();
    if let Some(cell) = cfd.resource_cell(TypeId::of::<ParcelSpec>()) {
        let mut b = cell.borrow_mut();
        let s = b.downcast_mut::<ParcelSpec>().unwrap();
        s.radius = radius;
        s.v_solid = v_solid;
    }
}

pub fn sph_bed_state(parent: &App) -> (f64, f64) {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let sph = subs.find("sph").unwrap();
    let atom_cell = sph.resource_cell(TypeId::of::<Atom>()).unwrap().borrow();
    let reg_cell = sph
        .resource_cell(TypeId::of::<AtomDataRegistry>())
        .unwrap()
        .borrow();
    let atoms = atom_cell.downcast_ref::<Atom>().unwrap();
    let registry = reg_cell.downcast_ref::<AtomDataRegistry>().unwrap();
    let sph = registry.expect::<SphAtom>("sph_bed_state");
    let n = atoms.nlocal as usize;
    let (mut sum_vz, mut sum_p) = (0.0, 0.0);
    let mut nfree = 0usize;
    for i in 0..n {
        if sph.is_boundary[i] > 0.5 {
            continue;
        }
        sum_vz += atoms.vel[i][2];
        sum_p += sph.pressure[i];
        nfree += 1;
    }
    let nf = nfree.max(1) as f64;
    (sum_vz / nf, sum_p / nf)
}

pub struct CoupledRunResult {
    pub mean_vz: f64,
    pub mean_p: f64,
    pub mom_err: f64,
}

pub fn coupled_diag(parent: &App) -> CoupledSeamDiag {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let cfd = subs.find("cfd").unwrap();
    let cell = cfd
        .resource_cell(TypeId::of::<CoupledSeamDiag>())
        .unwrap()
        .borrow();
    *cell.downcast_ref::<CoupledSeamDiag>().unwrap()
}

pub fn cleanup_subapps(parent: &mut App) {
    if let Some(cell) = parent.get_mut_resource(TypeId::of::<SubApps>()) {
        cell.borrow_mut()
            .downcast_mut::<SubApps>()
            .unwrap()
            .cleanup_all();
    }
}
