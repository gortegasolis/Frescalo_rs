//! Stage 2 (`frescalo`): local frequencies, effort standardisation,
//! benchmark species and sampling intensity.
//!
//! * **Local frequency** of a species at a site: the weighted share of the
//!   site's neighbourhood where it was recorded (any period).
//! * **Effort standardisation** ([`fresca`]): a multiplier α per site
//!   rescales the local frequencies, under the Poisson encounter model
//!   ([`crate::detection::poisson`]), until their frequency-weighted mean Φ
//!   equals the target (0.74 by default). Poorly recorded neighbourhoods need
//!   α > 1, well recorded ones α < 1.
//! * **Benchmark species**: at each site, the commonest species making up the
//!   top R\* of the expected species count. How many of them were recorded
//!   in a period measures the recording effort there and then (the
//!   *sampling intensity*).

use crate::config::{limits, FrescaloParams};
use crate::dataset::Dataset;
use crate::detection::poisson;
use crate::fortran::{binfnd, getnum, ld_i, ld_line, name_to_string, rec30_field, sort2, Arr2, Name};
use crate::output;
use std::io::Write;

/// Per-site sums over the weights file.
pub struct WeightTotals {
    /// Sum of neighbourhood weights for each target site.
    pub wgttot: Vec<f32>,
    /// Sum of squared weights for each target site.
    pub wgtt2: Vec<f32>,
    /// Number of targets each site is a neighbour of.
    pub numwgt: Vec<i32>,
}

/// Sum the neighbourhood weights of each site.
pub fn weight_totals(ds: &Dataset) -> WeightTotals {
    let mm = limits::FRESCALO_SITES;
    let mut t = WeightTotals {
        wgttot: vec![0.0f32; mm + 2],
        wgtt2: vec![0.0f32; mm + 2],
        numwgt: vec![0i32; mm + 2],
    };
    for iwgt in 1..=ds.nwgt {
        let dtji = ds.weights[iwgt];
        let f_neigh = rec30_field(&dtji, 0);
        let f_targ = rec30_field(&dtji, 1);
        let i = binfnd(&ds.sa, ds.m, &f_targ);
        let ii = binfnd(&ds.sa, ds.m, &f_neigh);
        let wgt = getnum(&rec30_field(&dtji, 2));
        t.wgttot[i] += wgt;
        t.wgtt2[i] += wgt * wgt;
        t.numwgt[ii] += 1;
    }
    t
}

/// Number of records (`idat`) and of species (`iitot`) at each site, from
/// the records sorted by (site, species, period).
pub fn records_per_site(ds: &Dataset) -> (Vec<i32>, Vec<i32>) {
    let mm = limits::FRESCALO_SITES;
    let mut idat = vec![0i32; mm + 2];
    let mut iitot = vec![0i32; mm + 2];
    let mut spec1 = crate::fortran::blank_name();
    for idtji in 1..=ds.nrec {
        let dtji = ds.records[idtji];
        let f_samp = rec30_field(&dtji, 0);
        let f_spec = rec30_field(&dtji, 1);
        let i = binfnd(&ds.sa, ds.m, &f_samp);
        idat[i] += 1;
        if f_spec != spec1 {
            iitot[i] += 1;
            spec1 = f_spec;
        }
    }
    (idat, iitot)
}

/// Local frequencies and presences.
pub struct LocalFrequencies {
    /// Local frequency of species `j` around site `i`.
    pub ffij: Arr2<f32>,
    /// 1 if species `j` was recorded at site `i` (any period).
    pub idata: Arr2<i32>,
    /// Value of the original's loop variable after the loop (`m + 1`), used
    /// by a progress-message quirk in the rescaling loop.
    pub ii_leftover: usize,
}

/// Pool each site's records over the neighbourhoods it belongs to.
pub fn local_frequencies(ds: &Dataset, idat: &[i32], totals: &WeightTotals) -> LocalFrequencies {
    let mm = limits::FRESCALO_SITES;
    let nn = limits::FRESCALO_SPECIES;
    let (m, n) = (ds.m, ds.n);
    let mut ffij: Arr2<f32> = Arr2::new(mm + 1, nn + 1, 0.0);
    let mut idata: Arr2<i32> = Arr2::new(mm + 1, nn + 1, 0);
    let mut jocc = vec![0i32; nn + 2];
    let mut idtji = 0usize;
    let mut iwgt = 0usize;
    let mut ii_leftover = 0usize;
    for ii in 1..=m {
        ii_leftover = ii + 1;
        if ii % 100 == 0 {
            ld_line(&format!("frequencies {}{}", name_to_string(&ds.sa[ii]), ld_i(ii as i64)));
        }
        for j in 1..=n {
            jocc[j] = 0;
            idata.set(ii, j, 0);
        }
        for _ in 1..=idat[ii] {
            idtji += 1;
            let f_spec = rec30_field(&ds.records[idtji], 1);
            let j = binfnd(&ds.sp, n, &f_spec);
            jocc[j] = 1;
            idata.set(ii, j, 1);
        }
        for _ in 1..=totals.numwgt[ii] {
            iwgt += 1;
            let dtji = ds.weights[iwgt];
            let i = binfnd(&ds.sa, m, &rec30_field(&dtji, 1));
            let wgt = getnum(&rec30_field(&dtji, 2));
            let denom = totals.wgttot[i] + 1.0E-10;
            for j in 1..=n {
                let contrib = (jocc[j] as f32) * wgt / denom;
                ffij.add(i, j, contrib);
            }
        }
    }
    LocalFrequencies { ffij, idata, ii_leftover }
}

/// Benchmark weights per species: 1, or 0.001 for species excluded from the
/// benchmarks. `exclusions` is 1-based.
pub fn benchmark_weights(ds: &Dataset, exclusions: &[Name]) -> Vec<f32> {
    let mut bwght = vec![1.0f32; limits::FRESCALO_SPECIES + 2];
    for name in exclusions.iter().skip(1) {
        let j = binfnd(&ds.sp, ds.n, name);
        if j != 0 {
            bwght[j] = 0.001;
        }
    }
    bwght
}

/// Result of standardising every site.
pub struct Standardised {
    /// Φ before rescaling (`Phi_in`) per site.
    pub phi_in: Vec<f32>,
    /// 1 if species `j` is a benchmark species at site `i`.
    pub ibench: Arr2<i32>,
    /// Total benchmark weight per site.
    pub abtot: Vec<f32>,
}

/// Run [`fresca`] for every site, replace the local frequencies in `lf.ffij`
/// by the rescaled ones, and pick each site's benchmark species.
#[allow(clippy::too_many_arguments)]
pub fn standardise_sites<W7: Write, W8: Write>(
    ds: &Dataset,
    lf: &mut LocalFrequencies,
    iitot: &[i32],
    totals: &WeightTotals,
    bwght: &[f32],
    params: &FrescaloParams,
    samples_out: &mut W7,
    frequencies_out: &mut W8,
) -> Standardised {
    let mm = limits::FRESCALO_SITES;
    let nn = limits::FRESCALO_SPECIES;
    let (m, n) = (ds.m, ds.n);
    let mut f = vec![0.0f32; nn + 2];
    let mut ff = vec![0.0f32; nn + 2];
    let mut jocc = vec![0i32; nn + 2];
    let mut jrank = vec![0i32; nn + 2];
    let mut out = Standardised {
        phi_in: vec![0.0f32; mm + 2],
        ibench: Arr2::new(mm + 1, nn + 1, 0),
        abtot: vec![0.0f32; mm + 2],
    };
    let ii = lf.ii_leftover;
    for i in 1..=m {
        // NB: the original mistakenly uses ii (leftover loop variable) here;
        // we reproduce the condition but avoid the out-of-bounds name lookup.
        if ii % 100 == 0 {
            let name = ds.sa.get(ii).map(name_to_string).unwrap_or_default();
            ld_line(&format!("rescaling {}{}", name, ld_i(ii as i64)));
        }
        for j in 1..=n {
            f[j] = lf.ffij.at(i, j);
            jocc[j] = lf.idata.at(i, j);
        }
        // wn2 is the effective number of weights in the neighbourhood
        let wn2 = totals.wgttot[i] * totals.wgttot[i] / (totals.wgtt2[i] + 1.0E-12);
        let (phi1, spnum) = fresca(
            i,
            n,
            iitot[i],
            &jocc,
            &mut f,
            &mut ff,
            &mut jrank,
            &ds.sa[i],
            &ds.sp,
            params,
            wn2,
            samples_out,
            frequencies_out,
        );
        out.phi_in[i] = phi1;
        out.abtot[i] = 1.0E-7;
        for j in 1..=n {
            if lf.ffij.at(i, j) != 0.0 {
                lf.ffij.set(i, j, ff[j]);
            }
            let jj = jrank[j] as usize;
            let rank1 = (j as f32) / spnum;
            // The case j=1 is included because with small samples rank1 may
            // be greater than the limit for j=1
            if rank1 < params.benchmark_limit || j == 1 {
                out.ibench.set(i, jj, 1);
                out.abtot[i] += bwght[jj];
            } else {
                out.ibench.set(i, jj, 0);
            }
        }
    }
    out
}

/// Sampling intensity per site and period, and the number of records per
/// species and period. Expects records ordered (species, period, site).
pub fn sampling_intensity(ds: &Dataset, st: &Standardised, bwght: &[f32]) -> (Arr2<f32>, Arr2<i32>) {
    let mut sampef: Arr2<f32> = Arr2::new(limits::FRESCALO_SITES + 1, limits::FRESCALO_PERIODS + 1, 1.0E-7);
    let mut lendat: Arr2<i32> = Arr2::new(limits::FRESCALO_SPECIES + 1, limits::FRESCALO_PERIODS + 1, 0);
    for idtji in 1..=ds.nrec {
        if idtji % 20000 == 0 {
            ld_line(&format!(
                "Main data to calc sampling effort{}{}",
                String::from_utf8_lossy(&ds.records[idtji]),
                ld_i(idtji as i64)
            ));
        }
        let dtji = ds.records[idtji];
        let j = binfnd(&ds.sp, ds.n, &rec30_field(&dtji, 0));
        let iit = binfnd(&ds.tim, ds.nt, &rec30_field(&dtji, 1));
        let i = binfnd(&ds.sa, ds.m, &rec30_field(&dtji, 2));
        lendat.add(j, iit, 1);
        let contrib = (st.ibench.at(i, j) as f32) * bwght[j] / st.abtot[i];
        sampef.add(i, iit, contrib);
    }
    (sampef, lendat)
}

/// The 98.5 percentile of the sites' input Φ, used to warn when the target Φ
/// looks too low. Sorts `phi_in[1..=m]` in place.
pub fn phi_in_985_percentile(phi_in: &mut [f32], m: usize) -> f32 {
    crate::fortran::sort_real(phi_in, m);
    let i985 = (0.985f32 * m as f32) as usize;
    phi_in[i985]
}

/// Effort standardisation for one site.
///
/// Finds the effort multiplier α such that, after each local frequency `f` is
/// rescaled to `1 − exp(−α·λ)` with `λ = −ln(1 − f)`, the frequency-weighted
/// mean frequency Φ = Σf'² / Σf' equals the target. α is reported in the
/// `Alpha` column of the samples file; a large α means the neighbourhood was
/// under-recorded. Also writes the site's samples-file row and its
/// frequencies-file rows, and leaves in `ff` the rescaled frequencies and in
/// `jrank` the species ordered by decreasing local frequency.
///
/// * `m` – site number, `n` – number of species, `itot` – species recorded at the site
/// * `jocc` – recorded (1) or not (0), per species
/// * `f` – local frequencies (clamped in place to `[fmin, fmax]`)
///
/// Returns (Φ before rescaling, expected species count after rescaling).
#[allow(clippy::too_many_arguments)]
pub fn fresca<W7: Write, W8: Write>(
    m: usize,
    n: usize,
    itot: i32,
    jocc: &[i32],
    f: &mut [f32],
    ff: &mut [f32],
    jrank: &mut [i32],
    samp1: &Name,
    splist: &[Name],
    params: &FrescaloParams,
    wn2: f32,
    unit7: &mut W7,
    unit8: &mut W8,
) -> (f32, f32) {
    let phibig = params.phi;
    let fmax = params.fmax;
    let fmin = params.fmin;
    let irepmx = params.max_iter;
    let mut alpha: f32 = 1.0;
    let mut phi: f32 = 0.0;
    let mut phi1: f32 = 0.0;
    let mut spnum: f32 = 0.0;
    let mut spnum1: f32 = 0.0;
    let mut ir = irepmx;
    let mut converged = false;
    for iter in 1..=irepmx {
        ir = iter;
        for j in 1..=n {
            if f[j] > fmax {
                f[j] = fmax;
            }
            if f[j] < fmin {
                f[j] = fmin;
            }
            ff[j] = poisson::intensity(f[j]);
        }
        let mut tot: f32 = 0.0;
        let mut tot2: f32 = 0.0;
        for j in 1..=n {
            // new frequency after recording intensity multiplied by alpha
            let ffij = poisson::prob_from_intensity(ff[j], alpha);
            tot += ffij;
            tot2 += ffij * ffij;
        }
        phi = tot2 / tot;
        spnum = tot;
        if iter < 20 {
            // successive approximation based on linear relation
            alpha = poisson::alpha_step(alpha, phi, phibig);
        } else {
            // crude successive approximation - slower with big datasets
            alpha = alpha * phibig / phi;
        }
        if iter == 1 {
            phi1 = phi;
            spnum1 = tot;
        }
        if (phi - phibig).abs() < params.tol {
            converged = true;
            break;
        }
    }
    // Fortran `DO ir=1,irepmx ... enddo` leaves ir = irepmx+1 after normal
    // (non-EXIT) completion; this only affects the reported Iter count on
    // non-convergence.
    if !converged {
        ir = irepmx + 1;
    }

    if m == 1 {
        output::writeln(unit7, output::SAMPLES_HEADER);
    }
    let alph = if alpha > 999.99 { 999.99 } else { alpha };
    output::write_sample_row(unit7, samp1, m, itot, phi1, alph, wn2, phi, spnum1, spnum, ir);

    for j in 1..=n {
        ff[j] = -f[j] + (j as f32) * 1.0E-12;
        jrank[j] = j as i32;
    }
    sort2(ff, jrank, n);
    for j in 1..=n {
        let jj = jrank[j] as usize;
        let mut fij = f[jj];
        if fij > fmax {
            fij = fmax;
        }
        let ffij = poisson::rescale(fij, alpha);
        let sdfij = (fij * (1.0 - fij) / wn2).sqrt();
        let mut ffff = fij + sdfij;
        if 1.0 - ffff < 1.0E-12 {
            ffff = 1.0 - 1.0E-12;
        }
        let fffff = fij - sdfij;
        let ffsd = poisson::rescale(ffff, alpha);
        let fffsd = poisson::rescale(fffff, alpha);
        // sdij is an estimate of the standard error
        let sdij = 0.5 * (ffsd - fffsd);
        if j == 1 && m == 1 {
            output::writeln(unit8, output::FREQUENCIES_HEADER);
        }
        ff[jj] = ffij;
        if f[jj] > 0.00005 {
            output::write_frequency_row(unit8, samp1, &splist[jj], jocc[jj], f[jj], ffij, sdij, j, (j as f32) / spnum);
        }
    }
    (phi1, spnum)
}
