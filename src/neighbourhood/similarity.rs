//! Stage 1b (`neighsim`, first half): reading the inputs and measuring how
//! alike sites are.
//!
//! Floristic similarity between two sites is Sørensen's coefficient on the
//! training-set species: twice the number of species they share, divided by
//! the sum of their species counts (0 = nothing in common, 1 = identical
//! lists). The training set is a well-recorded group used only to judge which
//! sites are ecologically alike.

use crate::config::limits;
use crate::fortran::{
    addwrd_guarded, binfnd, blank_name, isort, ld_i, make_rec30, name_to_string, rec30_field, sort30, Arr2, DataReader,
    Name, Rec30,
};

/// Everything `neighsim` reads: the `sampdist` ranking and the training set.
pub struct NeighsimInput {
    /// `sampdist` records (site, neighbour, spatial rank).
    pub distii: Vec<Rec30>,
    pub ndist: usize,
    /// Training-set records (species, site, unused), sorted by species.
    pub ddji: Vec<Rec30>,
    pub ndata: usize,
    /// Sites and species in the training set (sorted, 1-based).
    pub sa: Vec<Name>,
    pub m: usize,
    pub sp: Vec<Name>,
    pub n: usize,
    /// Species count per site and site count per species.
    pub itot: Vec<i32>,
    pub jtot: Vec<i32>,
    /// Last values of the reader's `samp`/`samp1` fields, reproduced in a
    /// later progress message exactly as the original prints them.
    pub last_samp: Name,
    pub last_samp1: Name,
}

/// Read the distance file, then the training set; sort the training records
/// and count species per site and sites per species. `Err` carries the
/// message for a fatal limit violation.
pub fn read_inputs(dist_reader: &mut DataReader, train_reader: &mut DataReader) -> Result<NeighsimInput, String> {
    let mm = limits::NEIGHSIM_SITES;
    let nn = limits::NEIGHSIM_SPECIES;
    let nndat = limits::NEIGHSIM_RECORDS;
    let mut sa = vec![blank_name(); mm + 2];
    let mut sp = vec![blank_name(); nn + 2];
    let mut ddji = vec![[b' '; 30]; nndat + 2];
    let mut distii = vec![[b' '; 30]; nndat + 2];
    let mut m: usize = 0;
    let mut n: usize = 0;
    let mut ndist: usize = 0;
    let mut ndata: usize = 0;

    // The same reader fields are shared by both loops, as in the original.
    let mut samp = blank_name();
    let mut samp1 = blank_name();
    let mut spec = blank_name();
    let mut any = blank_name();

    // Read in distance data
    while dist_reader.getd(&mut samp, &mut samp1, &mut any) {
        ndist += 1;
        if ndist > nndat {
            return Err(format!(
                " Too many data items in physical distance file - limit is{}",
                ld_i(nndat as i64)
            ));
        }
        if ndist % 20000 == 0 {
            crate::fortran::ld_line(&format!(
                "{} Dist {}{}{}",
                ld_i(ndist as i64),
                name_to_string(&samp),
                name_to_string(&samp1),
                name_to_string(&any)
            ));
        }
        distii[ndist] = make_rec30(&samp, &samp1, &any);
    }

    // Read in training-set species data
    while train_reader.getd(&mut samp, &mut spec, &mut any) {
        addwrd_guarded(&mut sa, mm, &mut m, &samp);
        if m > mm {
            return Err(format!(" Too many samples - limit is{}", ld_i(mm as i64)));
        }
        addwrd_guarded(&mut sp, nn, &mut n, &spec);
        if n > nn {
            return Err(format!(" Too many species - limit is{}", ld_i(nn as i64)));
        }
        ndata += 1;
        if ndata > nndat {
            return Err(format!(" Too many data items - limit is{}", ld_i(nndat as i64)));
        }
        if ndata % 20000 == 0 {
            crate::fortran::ld_line(&format!(
                "{} Spdata {}{}{}",
                ld_i(ndata as i64),
                name_to_string(&samp),
                name_to_string(&spec),
                name_to_string(&any)
            ));
        }
        ddji[ndata] = make_rec30(&spec, &samp, &any);
    }

    crate::fortran::cout(" Sorting main data ...");
    sort30(&mut ddji, ndata);
    crate::fortran::cout(" Sort completed");

    let mut itot = vec![0i32; mm + 2];
    let mut jtot = vec![0i32; nn + 2];
    for idata in 1..=ndata {
        if idata % 20000 == 0 {
            crate::fortran::ld_line(&format!(" Calculating totals{}", ld_i(idata as i64)));
        }
        let dji = ddji[idata];
        let i = binfnd(&sa, m, &rec30_field(&dji, 1));
        let j = binfnd(&sp, n, &rec30_field(&dji, 0));
        itot[i] += 1;
        jtot[j] += 1;
    }

    Ok(NeighsimInput {
        distii,
        ndist,
        ddji,
        ndata,
        sa,
        m,
        sp,
        n,
        itot,
        jtot,
        last_samp: samp,
        last_samp1: samp1,
    })
}

/// Number of training-set species shared by every pair of sites
/// (`simil(i1, i2)`; the diagonal holds each site's own species count).
/// `Err` carries the message if the sorted records are inconsistent.
pub fn shared_species(inp: &NeighsimInput) -> Result<Arr2<f32>, String> {
    let m = inp.m;
    let mm = limits::NEIGHSIM_SITES;
    let mut simil: Arr2<f32> = Arr2::new(mm + 1, mm + 1, 0.0);
    let mut iocc = vec![0i32; mm + 2];
    let mut idata = 0usize;
    for j in 1..=inp.n {
        for i in 1..=m {
            iocc[i] = 0;
        }
        for _ in 1..=inp.jtot[j] {
            idata += 1;
            if idata % 20000 == 0 {
                crate::fortran::ld_line(&format!(" Similarities{}", ld_i(idata as i64)));
            }
            let dji = inp.ddji[idata];
            let f_spec = rec30_field(&dji, 0);
            if f_spec != inp.sp[j] {
                return Err(format!(
                    "Unequal species{}{}",
                    name_to_string(&f_spec),
                    name_to_string(&inp.sp[j])
                ));
            }
            let i = binfnd(&inp.sa, m, &rec30_field(&dji, 1));
            iocc[i] = -(i as i32);
        }
        isort(&mut iocc, m);
        // mcc is the length of nonzero items in iocc
        let mut mcc = m;
        for ic in 1..=m {
            if iocc[ic] == 0 {
                mcc = ic - 1;
                break;
            }
        }
        for icc1 in 1..=mcc {
            let i1 = (-iocc[icc1]) as usize;
            for icc2 in 1..=mcc {
                let i2 = (-iocc[icc2]) as usize;
                simil.add(i1, i2, 1.0);
            }
        }
    }
    Ok(simil)
}

/// Sørensen similarity of site `i1` to every site, from the shared-species
/// counts: `2·shared / (species at i1 + species at i2)`, written into
/// `sim[1..=m]`.
pub fn sorensen_row(simil: &Arr2<f32>, itot: &[i32], i1: usize, m: usize, sim: &mut [f32]) {
    for i2 in 1..=m {
        sim[i2] = simil.at(i1, i2) * 2.0 / (itot[i1] + itot[i2]) as f32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorensen_on_toy_counts() {
        // Site 1 has 4 species, site 2 has 2, they share 2.
        let mut simil: Arr2<f32> = Arr2::new(2, 2, 0.0);
        simil.set(1, 1, 4.0);
        simil.set(2, 2, 2.0);
        simil.set(1, 2, 2.0);
        simil.set(2, 1, 2.0);
        let itot = [0, 4, 2];
        let mut sim = [0.0f32; 3];
        sorensen_row(&simil, &itot, 1, 2, &mut sim);
        assert_eq!(sim[1], 1.0);
        assert!((sim[2] - 2.0 * 2.0 / 6.0).abs() < 1e-7);
    }
}
