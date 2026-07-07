use cfd_ibm::coupling::drag_force_from_beta;
use field_core::{FvMesh, StructuredMesh, UniformMesh, UniformMeshConfig, Vec3};

use crate::drag::{beta_for, SeamMode};

/// One immersed parcel: center and solid volume `v_solid = m/rho_s`.
#[derive(Clone, Copy)]
pub struct Parcel {
    pub center: Vec3,
    pub v_solid: f64,
}

pub fn axis_centers(mesh: &UniformMesh) -> ([Vec<f64>; 3], usize) {
    let [ni, nj, nk] = mesh.dims();
    let ng = mesh.n_ghost();
    let xc = (0..ni)
        .map(|i| mesh.cell_centroid(mesh.idx_raw(i + ng, ng, ng))[0])
        .collect();
    let yc = (0..nj)
        .map(|j| mesh.cell_centroid(mesh.idx_raw(ng, j + ng, ng))[1])
        .collect();
    let zc = (0..nk)
        .map(|k| mesh.cell_centroid(mesh.idx_raw(ng, ng, k + ng))[2])
        .collect();
    ([xc, yc, zc], ng)
}

#[inline]
pub fn nearest_center(cs: &[f64], v: f64) -> usize {
    if cs.len() < 2 {
        return 0;
    }
    let dx = cs[1] - cs[0];
    (((v - cs[0]) / dx).round() as isize).clamp(0, cs.len() as isize - 1) as usize
}

pub fn containing_cell(mesh: &UniformMesh, centers: &[Vec<f64>; 3], ng: usize, p: Vec3) -> usize {
    let i = nearest_center(&centers[0], p[0]);
    let j = nearest_center(&centers[1], p[1]);
    let k = nearest_center(&centers[2], p[2]);
    mesh.idx_raw(i + ng, j + ng, k + ng)
}

/// Per-cell void fraction `eps = 1 - sum(V_solid)/V_cell` by containment
/// deposition, plus each parcel's containing-cell index.
pub fn deposit_bed_void_fraction(mesh: &UniformMesh, parcels: &[Parcel]) -> (Vec<f64>, Vec<usize>) {
    let (centers, ng) = axis_centers(mesh);
    let total = mesh.n_cells_total();
    let mut solid = vec![0.0f64; total];
    let mut cell_of_parcel = Vec::with_capacity(parcels.len());
    for p in parcels {
        let c = containing_cell(mesh, &centers, ng, p.center);
        solid[c] += p.v_solid;
        cell_of_parcel.push(c);
    }
    let mut eps = vec![1.0f64; total];
    for c in 0..total {
        let v = mesh.cell_volume(c);
        if v > 0.0 {
            eps[c] = (1.0 - solid[c] / v).clamp(1e-6, 1.0);
        }
    }
    (eps, cell_of_parcel)
}

#[allow(clippy::too_many_arguments)]
pub fn build_deposit_mesh(
    x_lo: f64,
    x_hi: f64,
    y_lo: f64,
    y_hi: f64,
    z_bottom: f64,
    z_height: f64,
    cell: f64,
) -> UniformMeshConfig {
    let nx = (((x_hi - x_lo) / cell).round() as usize).max(1);
    let ny = (((y_hi - y_lo) / cell).round() as usize).max(1);
    UniformMeshConfig {
        nx: nx + 2,
        ny: ny + 2,
        nz: 1,
        ng: 1,
        bounds_lo: [x_lo - cell, y_lo - cell, z_bottom],
        bounds_hi: [
            x_lo + (nx + 1) as f64 * cell,
            y_lo + (ny + 1) as f64 * cell,
            z_bottom + z_height,
        ],
        y_edges: None,
        z_edges: None,
    }
}

#[allow(clippy::too_many_arguments)]
pub fn bed_fluid_force_z(
    parcels: &[Parcel],
    eps_field: &[f64],
    cell_of_parcel: &[usize],
    eps_bed: f64,
    u_super: f64,
    rho_f: f64,
    mu: f64,
    d_grain: f64,
    g: f64,
    mode: SeamMode,
) -> f64 {
    let u_g = u_super / eps_bed;
    let rel: Vec3 = [0.0, 0.0, u_g];
    let rel_speed = u_g;
    let mut f_fluid = 0.0;
    for (i, p) in parcels.iter().enumerate() {
        let eps = eps_field[cell_of_parcel[i]];
        let beta = beta_for(mode, eps, rho_f, mu, d_grain, rel_speed);
        let drag = drag_force_from_beta(beta, p.v_solid, eps, rel);
        let pg_coeff = if mode.omit_pressure_grad {
            0.0
        } else {
            p.v_solid * beta / eps
        };
        let buoy = rho_f * p.v_solid * g;
        let mut fz = drag[2] + pg_coeff * rel[2] + buoy;
        if mode.corrupt_eps_power {
            fz /= eps_bed;
        }
        f_fluid += fz;
    }
    f_fluid
}

#[allow(clippy::too_many_arguments)]
pub fn measure_umf(
    parcels: &[Parcel],
    eps_field: &[f64],
    cell_of_parcel: &[usize],
    eps_bed: f64,
    rho_f: f64,
    rho_s: f64,
    mu: f64,
    d_grain: f64,
    g: f64,
    mode: SeamMode,
) -> f64 {
    let w_full: f64 = parcels.iter().map(|p| rho_s * p.v_solid * g).sum();
    let net = |u: f64| {
        bed_fluid_force_z(
            parcels,
            eps_field,
            cell_of_parcel,
            eps_bed,
            u,
            rho_f,
            mu,
            d_grain,
            g,
            mode,
        ) - w_full
    };
    let (mut lo, mut hi) = (1e-5, 50.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if net(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn deposit_worst_err(mesh: &UniformMesh, eps_field: &[f64], eps_bed: f64) -> f64 {
    let mut e = 0.0f64;
    for (c, &eps) in eps_field.iter().enumerate() {
        if !mesh.is_local_cell(c) || eps >= 1.0 - 1e-9 {
            continue;
        }
        e = e.max((eps - eps_bed).abs() / eps_bed);
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn force_balance_measurement_uses_drag_pressure_gradient_and_buoyancy() {
        let parcels = vec![Parcel {
            center: [0.0, 0.0, 0.0],
            v_solid: 2.0e-9,
        }];
        let eps_field = vec![0.4];
        let cells = vec![0usize];
        let (eps, rho_f, rho_s, mu, d, g) = (0.4, 1.2, 2500.0, 1.8e-5, 5.0e-4, 9.81);
        let u = measure_umf(
            &parcels,
            &eps_field,
            &cells,
            eps,
            rho_f,
            rho_s,
            mu,
            d,
            g,
            SeamMode::default(),
        );
        let net = bed_fluid_force_z(
            &parcels,
            &eps_field,
            &cells,
            eps,
            u,
            rho_f,
            mu,
            d,
            g,
            SeamMode::default(),
        ) - rho_s * parcels[0].v_solid * g;
        assert!(net.abs() / (rho_s * parcels[0].v_solid * g) < 1e-12);
    }
}
