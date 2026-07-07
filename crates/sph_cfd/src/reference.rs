pub fn archimedes(rho_f: f64, rho_s: f64, g: f64, d: f64, mu: f64) -> f64 {
    rho_f * (rho_s - rho_f) * g.abs() * d.powi(3) / (mu * mu)
}

/// Wen & Yu (1966): `Re_mf = sqrt(33.7² + 0.0408 Ar) - 33.7`,
/// `U_mf = Re_mf μ/(ρ_f d)`.
pub fn u_mf_wen_yu(rho_f: f64, rho_s: f64, g: f64, d: f64, mu: f64) -> (f64, f64, f64) {
    let ar = archimedes(rho_f, rho_s, g, d, mu);
    let re_mf = (33.7f64 * 33.7 + 0.0408 * ar).sqrt() - 33.7;
    (re_mf * mu / (rho_f * d), ar, re_mf)
}

/// Superficial velocity where a packed-bed pressure drop (`c1` viscous, `c2`
/// inertial) equals the buoyant weight per length `(1-eps)(rho_s-rho_f)g`.
#[allow(clippy::too_many_arguments)]
pub fn u_mf_balance(
    c1: f64,
    c2: f64,
    eps: f64,
    rho_f: f64,
    rho_s: f64,
    g: f64,
    d: f64,
    mu: f64,
) -> f64 {
    let om = 1.0 - eps;
    let e3 = eps.powi(3);
    let a_visc = c1 * om / e3 * mu / (d * d);
    let a_inert = c2 / e3 * rho_f / d;
    let target = (rho_s - rho_f) * g.abs();
    (-a_visc + (a_visc * a_visc + 4.0 * a_inert * target).sqrt()) / (2.0 * a_inert)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn force_balance_solution_satisfies_packed_bed_quadratic() {
        let (c1, c2, eps, rho_f, rho_s, g, d, mu) =
            (180.0, 1.8, 0.41, 1.2, 2650.0, 9.81, 5.0e-4, 1.8e-5);
        let u = u_mf_balance(c1, c2, eps, rho_f, rho_s, g, d, mu);
        let om = 1.0 - eps;
        let lhs = c1 * om / eps.powi(3) * mu * u / (d * d) + c2 / eps.powi(3) * rho_f * u * u / d;
        let rhs = (rho_s - rho_f) * g;
        assert!((lhs - rhs).abs() / rhs < 1e-12);
    }
}
