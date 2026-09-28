//! Hill's time factor for one species in one period (`tfcalc`).
//!
//! The expected chance of recording the species at site *i* in period *t* is
//! `P_it = 1 − exp(−λ_i · x_t)`, where `λ_i` is the encounter intensity from
//! the standardised local frequency times the sampling intensity, and `x_t`
//! is the time factor. Starting from `x_t = 1`, the time factor is rescaled
//! until the expected number of sites with records equals the observed
//! number. Its standard error is found by repeating the search with the
//! observed number raised by one standard deviation.
//!
//! Sites with a sampling intensity below 0.0995 (no systematic recording)
//! are down-weighted.

use crate::detection::poisson;

/// Time factor and its supporting counts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimeFactor {
    /// Time factor `x_t`.
    pub tf: f32,
    /// Standard error of `tf`.
    pub sd: f32,
    /// Weighted observed number of sites with records.
    pub sptot: f32,
    /// Unweighted number of sites with records.
    pub jtot: i32,
    /// Weighted expected number of sites with records.
    pub esttot: f32,
    /// Sites with a non-zero chance of a record.
    pub ic1: i32,
    /// Sites where that chance was capped at 0.98.
    pub ic2: i32,
}

/// Compute the time factor.
///
/// * `iocc[i]` – 1 if the species was recorded at site `i` in the period
/// * `smpint[i]` – sampling intensity at site `i` in the period
/// * `fff[i]` – standardised local frequency of the species at site `i`
pub fn tfcalc(iocc: &[i32], smpint: &[f32], fff: &[f32], m: usize) -> TimeFactor {
    const KMAX: i32 = 100;
    let mut tf: f32 = 1.0;
    let mut esttot: f32 = 0.0;
    let mut estvar: f32 = 0.0;
    let mut sptot: f32 = 0.0;
    let mut jtot: i32 = 0;
    let mut ic1: i32 = 0;
    let mut ic2: i32 = 0;
    for _ in 1..=KMAX {
        esttot = 0.0;
        estvar = 0.0;
        sptot = 0.0;
        jtot = 0;
        ic1 = 0;
        ic2 = 0;
        for i in 1..=m {
            let mut wgt: f32 = 1.0;
            if smpint[i] < 0.0995 {
                wgt = 10.0 * smpint[i] + 0.005;
            }
            // probability of finding the species, times the sampling intensity
            let mut pfac = smpint[i] * fff[i];
            if pfac > 0.0 {
                ic1 += 1;
            }
            if pfac > 0.98 {
                // cap so that the log can be taken
                pfac = 0.98;
                ic2 += 1;
            }
            let plog = poisson::intensity(pfac);
            let estval = poisson::prob_from_intensity(plog, tf);
            esttot += wgt * estval;
            estvar += wgt * wgt * estval * (1.0 - estval);
            sptot += wgt * (iocc[i] as f32);
            jtot += iocc[i];
        }
        if (sptot - esttot).abs() < 0.0005 {
            break;
        }
        tf = tf * sptot / (esttot + 0.0000001);
    }

    // sptot1 is precisely 1 standard deviation bigger than sptot;
    // recalculate with this target value to obtain the standard error of tf.
    let sptot1 = sptot + estvar.sqrt();
    let mut tf1 = tf;
    for _ in 1..=KMAX {
        let mut esttt1: f32 = 0.0;
        for i in 1..=m {
            let mut wgt: f32 = 1.0;
            if smpint[i] < 0.0995 {
                wgt = 10.0 * smpint[i] + 0.005;
            }
            let mut pfac = smpint[i] * fff[i];
            if pfac > 0.98 {
                pfac = 0.98;
            }
            let plog = poisson::intensity(pfac);
            let estval = poisson::prob_from_intensity(plog, tf1);
            esttt1 += wgt * estval;
        }
        if (sptot1 - esttt1).abs() < 0.0005 {
            break;
        }
        tf1 = tf1 * sptot1 / (esttt1 + 0.0000001);
    }

    TimeFactor {
        tf,
        sd: tf1 - tf,
        sptot,
        jtot,
        esttot,
        ic1,
        ic2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_expectation_gives_unit_time_factor() {
        // Two sites, fully sampled, frequency 0.5 each, one recorded:
        // expected = 2 × 0.5 = 1 = observed, so tf stays 1.
        let iocc = [0, 1, 0];
        let smpint = [0.0, 1.0, 1.0];
        let fff = [0.0, 0.5, 0.5];
        let r = tfcalc(&iocc, &smpint, &fff, 2);
        assert!((r.tf - 1.0).abs() < 1e-3);
        assert_eq!(r.jtot, 1);
        assert_eq!(r.ic1, 2);
        assert_eq!(r.ic2, 0);
        assert!(r.sd > 0.0);
    }

    #[test]
    fn more_records_than_expected_gives_increase() {
        let iocc = [0, 1, 1, 1, 0];
        let smpint = [0.0, 1.0, 1.0, 1.0, 1.0];
        let fff = [0.0, 0.4, 0.4, 0.4, 0.4];
        assert!(tfcalc(&iocc, &smpint, &fff, 4).tf > 1.0);
    }
}
