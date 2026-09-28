//! Stage 1b (`neighsim`, second half): choosing each site's neighbourhood
//! and weighting its members.
//!
//! Only sites on the spatial shortlist from `sampdist` are eligible. Among
//! them, the *K* most floristically similar form the neighbourhood (Hill 2012
//! used K = 100). Each member gets the weight
//!
//! `w = (1 − r_sim²)⁴ × (1 − r_dist²)⁴`
//!
//! where `r_sim` and `r_dist` are its similarity rank and distance rank scaled
//! to run from 0 (the site itself) towards 1 (the edge of the neighbourhood).
//! The site itself gets weight 1; weights fall off smoothly with rank rather
//! than stopping abruptly, and weights of 0.00005 or less are dropped.

use super::similarity::{sorensen_row, NeighsimInput};
use crate::config::limits;
use crate::fortran::{binfnd, getnum, ld_i, ld_line, name_to_string, rec30_field, sort2, Arr2};
use crate::output;
use std::io::Write;

/// Factor that lifts every site on the spatial shortlist above every site
/// that is not, so that similarity ranks only ever pick shortlisted sites.
pub const BIG: f32 = 1_000_000.0;

/// Weights at or below this are not written.
pub const SMALL: f32 = 0.00005;

/// The smooth, bounded rank kernel `(1 − t²)⁴`: 1 at t = 0, 0 at t = 1.
#[inline]
pub fn rank_kernel(t: f32) -> f32 {
    (1.0 - t * t).powi(4)
}

/// Mark the spatial shortlist: multiply the shared-species count of every
/// (site, shortlisted neighbour) pair by [`BIG`] and record the neighbour's
/// spatial rank in `iseqq`. Returns the largest spatial rank seen (the
/// shortlist size *N* used by `sampdist`), or 0 if no site names matched.
pub fn mark_spatial_shortlist(inp: &NeighsimInput, simil: &mut Arr2<f32>) -> (Arr2<i32>, i32) {
    let mm = limits::NEIGHSIM_SITES;
    let mut iseqq: Arr2<i32> = Arr2::new(mm + 1, mm + 1, 0);
    let mut neigh1: i32 = 0;
    let mut samp = inp.last_samp;
    let mut samp1 = inp.last_samp1;
    for idist in 1..=inp.ndist {
        if idist % 20000 == 0 {
            // NB: the original prints samp/samp1 from the *previous* record here.
            ld_line(&format!(
                "{} Dist {}{}",
                ld_i(idist as i64),
                name_to_string(&samp),
                name_to_string(&samp1)
            ));
        }
        let dji = inp.distii[idist];
        samp = rec30_field(&dji, 0);
        samp1 = rec30_field(&dji, 1);
        let i1 = binfnd(&inp.sa, inp.m, &samp);
        let i2 = binfnd(&inp.sa, inp.m, &samp1);
        if i1 != 0 && i2 != 0 {
            simil.set(i1, i2, simil.at(i1, i2) * BIG);
            let seq = getnum(&rec30_field(&dji, 2));
            let iseq = seq as i32; // ifix: truncation towards zero
            iseqq.set(i1, i2, iseq);
            if neigh1 < iseq {
                neigh1 = iseq;
            }
        }
    }
    (iseqq, neigh1)
}

/// For every site: rank all sites by similarity (shortlisted first), write
/// the top `neigh` to the similarity file and, with their combined weights,
/// to the weights file.
pub fn write_neighbourhoods<W1: Write, W2: Write>(
    inp: &NeighsimInput,
    simil: &Arr2<f32>,
    iseqq: &Arr2<i32>,
    neigh: i32,
    neigh1: i32,
    sim_out: &mut W1,
    wgt_out: &mut W2,
) {
    let m = inp.m;
    let mut sim = vec![0.0f32; m + 2];
    let mut index = vec![0i32; m + 2];
    for i1 in 1..=m {
        if i1 % 100 == 0 {
            ld_line(&format!("Writing output  {}{}", name_to_string(&inp.sa[i1]), ld_i(i1 as i64)));
        }
        sorensen_row(simil, &inp.itot, i1, m, &mut sim);
        for (i2, ix) in index.iter_mut().enumerate().take(m + 1).skip(1) {
            *ix = i2 as i32;
        }
        sort2(&mut sim, &mut index, m);
        for is2 in 1..=(neigh.max(0) as usize) {
            // most similar first: walk the ascending sort from the top
            let i2 = m as isize - is2 as isize + 1;
            if i2 < 1 {
                break;
            }
            let i2 = i2 as usize;
            let iis2 = index[i2] as usize;
            output::write_similarity_row(sim_out, &inp.sa[i1], &inp.sa[iis2], is2 as i64, sim[i2] / BIG);
            if neigh1 == 0 {
                continue;
            }
            let t = (is2 as f32 - 1.0) / neigh as f32;
            let amult1 = rank_kernel(t);
            let t2 = (iseqq.at(i1, iis2) - 1) as f32 / neigh1 as f32;
            let amult2 = rank_kernel(t2);
            let amult = amult1 * amult2;
            if amult > SMALL {
                output::write_weight_row(wgt_out, &inp.sa[i1], &inp.sa[iis2], amult, amult1, amult2, neigh, neigh1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kernel_endpoints_and_shape() {
        assert_eq!(rank_kernel(0.0), 1.0);
        assert_eq!(rank_kernel(1.0), 0.0);
        assert!(rank_kernel(0.25) > rank_kernel(0.5));
        assert!(rank_kernel(0.5) > rank_kernel(0.75));
    }
}
