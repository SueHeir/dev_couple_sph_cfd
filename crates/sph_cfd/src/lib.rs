//! # sph_cfd — reusable granular-SPH ↔ CFD coupling seam helpers
//!
//! This crate holds the general pieces every unresolved granular-SPH/CFD case repeats:
//! declarative gas/grid/gravity blocks, packed-bed closure/reference math, parcel
//! deposition and force-balance helpers, and the CFD-side `grass_multi` seam systems.
//! Case geometry, validation tolerances, comparison cases, and plots stay in `examples/`.

pub mod bed;
pub mod config;
pub mod drag;
pub mod reference;
#[cfg(feature = "mpi-routing")]
pub mod routing;
pub mod seam;

pub mod prelude {
    pub use crate::bed::{
        axis_centers, bed_fluid_force_z, build_deposit_mesh, containing_cell,
        deposit_bed_void_fraction, deposit_worst_err, measure_umf, nearest_center, Parcel,
    };
    pub use crate::config::{GasCfg, GravityCfg, GridCfg};
    pub use crate::drag::{beta_for, ergun_beta, macdonald_beta, SeamMode};
    pub use crate::reference::{archimedes, u_mf_balance, u_mf_wen_yu};
    pub use crate::seam::{
        add_standard_coupled_schedule, build_coupled_cfd, cleanup_subapps, coupled_diag,
        coupled_seam_system, export_kinematics, import_forces, prime_sph_and_spec, read_sph_bed,
        sph_bed_state, sph_fluid_force, CoupledGasProps, CoupledRunResult, CoupledSeamDiag,
        FluidForces, ParcelSpec, Phase, R_GAS,
    };
}
