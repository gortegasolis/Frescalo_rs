//! The Poisson encounter model.
//!
//! Picture a recorder searching a site. Suppose encounters with a species
//! happen at random, at a steady average rate, so that the number of
//! encounters during a visit follows a Poisson distribution with mean λ (the
//! *intensity*: how often the species would be met with the effort spent).
//! The species gets recorded if it is met at least once, so
//!
//! `P(recorded) = 1 − e^(−λ)`.
//!
//! Doubling the effort doubles λ. In general, multiplying effort by `a` turns
//! a recording probability `p` into
//!
//! `p' = 1 − e^(−a·λ) = 1 − (1 − p)^a`,  with `λ = −ln(1 − p)`.
//!
//! Two parts of FRESCALO rely on this:
//!
//! * **Effort standardisation** (`fresca`): each neighbourhood's local
//!   frequencies are rescaled with a multiplier α so that every site is
//!   compared at the same standard effort.
//! * **Time factors** (`tfcalc`): the chance of recording a species at a site
//!   in a period is its intensity there multiplied by the time factor.
//!
//! Statistically this is a binomial model with a complementary log-log link:
//! `ln(−ln(1 − p))` shifts by `ln a` when effort is multiplied by `a`.
//!
//! **Why two functions compute the same thing.** [`prob_from_intensity`] and
//! [`rescale`] are mathematically identical, but the original code evaluates
//! them in different orders of single-precision operations, which can differ
//! in the last digit. Both forms are kept, each used exactly where the
//! original used it, so that results stay byte-identical with the Fortran.

/// Encounter intensity λ from a recording probability: `−ln(1 − p)`.
#[inline]
pub fn intensity(p: f32) -> f32 {
    -(1.0 - p).ln()
}

/// Recording probability when intensity `lambda` is multiplied by `a`,
/// evaluated as `1 − exp(−λ·a)` (used in the α search and in `tfcalc`).
#[inline]
pub fn prob_from_intensity(lambda: f32, a: f32) -> f32 {
    1.0 - (-lambda * a).exp()
}

/// Recording probability `p` after multiplying effort by `a`, evaluated as
/// `1 − exp(a·ln(1 − p))` (used for the rescaled frequencies and their
/// standard errors).
#[inline]
pub fn rescale(p: f32, a: f32) -> f32 {
    1.0 - (a * (1.0 - p).ln()).exp()
}

/// Hill's accelerated update of the effort multiplier α, used in the first
/// 19 iterations of the α search. It exploits the near-linear relation
/// between `ln(1 − Φ)` and `ln α` under this model; 1.86 is Hill's empirical
/// slope, specific to the Poisson model.
#[inline]
pub fn alpha_step(alpha: f32, phi: f32, target: f32) -> f32 {
    alpha * (1.86f32 * ((1.0 - phi).ln() - (1.0 - target).ln())).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_effort_leaves_probability_unchanged() {
        for &p in &[0.01f32, 0.2, 0.5, 0.74, 0.99] {
            assert!((rescale(p, 1.0) - p).abs() < 1e-6);
            assert!((prob_from_intensity(intensity(p), 1.0) - p).abs() < 1e-6);
        }
    }

    #[test]
    fn more_effort_means_higher_probability() {
        let p = 0.3f32;
        assert!(rescale(p, 0.5) < p);
        assert!(rescale(p, 2.0) > p);
        assert!(rescale(p, 4.0) > rescale(p, 2.0));
    }

    #[test]
    fn two_evaluation_forms_agree_closely() {
        for &p in &[0.001f32, 0.1, 0.5, 0.9, 0.999] {
            for &a in &[0.2f32, 1.0, 3.7, 20.0] {
                let x = prob_from_intensity(intensity(p), a);
                let y = rescale(p, a);
                assert!((x - y).abs() <= 4.0 * f32::EPSILON, "p={} a={} {} {}", p, a, x, y);
            }
        }
    }

    #[test]
    fn alpha_step_moves_towards_target() {
        // Φ below target → α grows; Φ above target → α shrinks; at target → unchanged.
        assert!(alpha_step(1.0, 0.6, 0.74) > 1.0);
        assert!(alpha_step(1.0, 0.8, 0.74) < 1.0);
        assert_eq!(alpha_step(1.5, 0.74, 0.74), 1.5);
    }
}
