/// MacDonald et al. (1979) interphase coefficient β (the "Ergun revisited" 180/1.8
/// re-fit).
pub fn macdonald_beta(eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    let eps = eps.clamp(1e-6, 1.0);
    let om = 1.0 - eps;
    180.0 * om * om * mu / (eps * d * d) + 1.8 * om * rho_f * rel_speed / d
}

/// Ergun (1952) β (150/1.75).
pub fn ergun_beta(eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    let eps = eps.clamp(1e-6, 1.0);
    let om = 1.0 - eps;
    150.0 * om * om * mu / (eps * d * d) + 1.75 * om * rho_f * rel_speed / d
}

/// Which β closure the seam assembles, and whether to inject a fault for negative
/// controls.
#[derive(Clone, Copy)]
pub struct SeamMode {
    /// `true` -> MacDonald(1979) measured closure; `false` -> Ergun(1952).
    pub macdonald: bool,
    /// Drop the ∇P pressure-gradient force.
    pub omit_pressure_grad: bool,
    /// Inject the ε²-instead-of-ε³ force-reduction bug.
    pub corrupt_eps_power: bool,
}

impl Default for SeamMode {
    fn default() -> Self {
        Self {
            macdonald: true,
            omit_pressure_grad: false,
            corrupt_eps_power: false,
        }
    }
}

pub fn beta_for(mode: SeamMode, eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    if mode.macdonald {
        macdonald_beta(eps, rho_f, mu, d, rel_speed)
    } else {
        ergun_beta(eps, rho_f, mu, d, rel_speed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn beta_closures_match_published_coefficients() {
        let (eps, rho, mu, d, u) = (0.42, 1.2, 1.8e-5, 4.0e-4, 0.73);
        let om = 1.0 - eps;
        let mac = 180.0 * om * om * mu / (eps * d * d) + 1.8 * om * rho * u / d;
        let erg = 150.0 * om * om * mu / (eps * d * d) + 1.75 * om * rho * u / d;
        assert!((macdonald_beta(eps, rho, mu, d, u) - mac).abs() < 1e-12 * mac);
        assert!((ergun_beta(eps, rho, mu, d, u) - erg).abs() < 1e-12 * erg);
        assert!((beta_for(SeamMode::default(), eps, rho, mu, d, u) - mac).abs() < 1e-12 * mac);
    }
}
