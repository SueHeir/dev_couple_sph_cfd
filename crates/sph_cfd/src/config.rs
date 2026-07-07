use serde::Deserialize;

#[derive(Deserialize, Default)]
pub struct GasCfg {
    pub rho: f64,
    pub p: f64,
    pub mu: f64,
}

#[derive(Deserialize, Default)]
pub struct GridCfg {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub ng: usize,
    pub z_hi: f64,
}

#[derive(Deserialize, Default)]
pub struct GravityCfg {
    pub gz: f64,
}
