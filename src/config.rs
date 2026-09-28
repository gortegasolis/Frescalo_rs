//! Parameters of the method and their defaults, plus the array limits
//! inherited from the Fortran originals.
//!
//! The defaults are the values Hill (2012) used for British data: a spatial
//! shortlist of 200 nearest sites, neighbourhoods of the 100 most similar of
//! those, a target local frequency Φ = 0.74 and a benchmark limit R\* = 0.2703.

use crate::distance::Metric;

/// Default number of nearest sites `sampdist` keeps per site (the spatial
/// shortlist, *N*).
pub const DEFAULT_SPATIAL_NEIGHBOURS: i32 = 200;

/// Default number of neighbours `neighsim` keeps per site, chosen by
/// floristic similarity from the spatial shortlist (*K*).
pub const DEFAULT_NEIGHBOURHOOD_SIZE: i32 = 100;

/// Settings for `sampdist`.
#[derive(Clone, Copy, Debug)]
pub struct SampdistParams {
    /// Number of nearest sites written per site (including the site itself).
    pub neighbours: i32,
    /// How distances between sites are measured.
    pub metric: Metric,
}

impl Default for SampdistParams {
    fn default() -> Self {
        SampdistParams {
            neighbours: DEFAULT_SPATIAL_NEIGHBOURS,
            metric: Metric::Planar,
        }
    }
}

/// Settings for `neighsim`.
#[derive(Clone, Copy, Debug)]
pub struct NeighsimParams {
    /// Number of neighbours kept per site, ranked by floristic similarity.
    pub neighbours: i32,
}

impl Default for NeighsimParams {
    fn default() -> Self {
        NeighsimParams {
            neighbours: DEFAULT_NEIGHBOURHOOD_SIZE,
        }
    }
}

/// Settings for `frescalo`.
#[derive(Clone, Copy, Debug)]
pub struct FrescaloParams {
    /// Target local frequency Φ that every neighbourhood is standardised to.
    pub phi: f32,
    /// Benchmark limit R\*: species ranked in the top R\* of a site's expected
    /// species count are benchmark species.
    pub benchmark_limit: f32,
    /// Largest local frequency allowed (< 1 so that ln(1 − f) exists).
    pub fmax: f32,
    /// Smallest local frequency allowed (keeps frequencies from summing to zero).
    pub fmin: f32,
    /// Convergence tolerance on Φ when searching for the effort multiplier α.
    pub tol: f32,
    /// Maximum number of iterations in the α search.
    pub max_iter: i32,
}

/// Default target Φ.
pub const DEFAULT_PHI: f32 = 0.74;
/// Default benchmark limit R\*.
pub const DEFAULT_BENCHMARK_LIMIT: f32 = 0.2703;
/// Accepted range for Φ.
pub const PHI_RANGE: (f32, f32) = (0.50, 0.95);
/// Accepted range for R\*.
pub const BENCHMARK_LIMIT_RANGE: (f32, f32) = (0.08, 0.5);

impl Default for FrescaloParams {
    fn default() -> Self {
        FrescaloParams {
            phi: DEFAULT_PHI,
            benchmark_limit: DEFAULT_BENCHMARK_LIMIT,
            fmax: 0.99999,
            fmin: 1.0E-10,
            tol: 0.0003,
            max_iter: 100,
        }
    }
}

/// Array limits of the original programs (Fortran `parameter` values).
pub mod limits {
    /// `sampdist`: maximum number of locations.
    pub const SAMPDIST_SITES: usize = 400_000;

    /// `neighsim`: maximum number of sites.
    pub const NEIGHSIM_SITES: usize = 4000;
    /// `neighsim`: maximum number of training-set species.
    pub const NEIGHSIM_SPECIES: usize = 10_000;
    /// `neighsim`: maximum number of training-set records (and distance records).
    pub const NEIGHSIM_RECORDS: usize = 5_000_000;

    /// `frescalo`: maximum number of sites.
    pub const FRESCALO_SITES: usize = 4000;
    /// `frescalo`: maximum number of species.
    pub const FRESCALO_SPECIES: usize = 2000;
    /// `frescalo`: maximum number of time periods.
    pub const FRESCALO_PERIODS: usize = 100;
    /// `frescalo`: maximum number of occurrence records.
    pub const FRESCALO_RECORDS: usize = 2_000_000;
    /// `frescalo`: maximum number of neighbourhood weights.
    pub const FRESCALO_WEIGHTS: usize = 500_000;
}
