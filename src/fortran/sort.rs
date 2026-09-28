//! Sorting and sorted-list lookup (Fortran heapsorts, `binfnd`, `addwrd`).

use super::text::{Name, Rec30};

// ---------------------------------------------------------------------------
// binfnd: binary search in a sorted 1-based list of names.
// Returns the 1-based index, or 0 if not found.
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself. Looks up a site or species code in a sorted list and
/// returns its position (or 0 if not found), using binary search so this stays fast even with
/// millions of records — e.g. the 2 038 000 bryophyte records across 3695 sites used in the
/// paper's example (Hill 2012, p.196). Every time the code needs to relate a record to "which site
/// is this" or "which species is this", it goes through `binfnd`.
pub fn binfnd(ma: &[Name], n: usize, na: &Name) -> usize {
    let mut imin = 1usize;
    let mut iamin = ma[imin];
    let mut imax = n;
    let mut iamax = ma[imax];
    loop {
        if imax - imin <= 1 {
            if iamin == *na {
                return imin;
            }
            if iamax == *na {
                return imax;
            }
            return 0;
        }
        let imid = (imax + imin) / 2;
        let iamid = ma[imid];
        if *na <= iamid {
            imax = imid;
            iamax = iamid;
        } else {
            imin = imid;
            iamin = iamid;
        }
    }
}

// ---------------------------------------------------------------------------
// addwrd: append samp to the sorted list sa (1-based, length *m) if new.
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself. Builds up the master, sorted list of distinct site or
/// species codes seen so far as records are read in one at a time — this is how the programs learn
/// the full set of sites/species in a data set without being told it in advance. `addwrd_guarded`
/// is the same thing but stops (silently) once a fixed array-size limit is reached, matching the
/// original Fortran's static array bounds.
pub fn addwrd(sa: &mut [Name], m: &mut usize, samp: &Name) {
    let i = if *m == 0 { 0 } else { binfnd(sa, *m, samp) };
    if i != 0 {
        return;
    }
    *m += 1;
    sa[*m] = *samp;
    sa[1..=*m].sort_unstable();
}

/// neighsim variant of addwrd with the extra overflow guard.
pub fn addwrd_guarded(sa: &mut [Name], mm: usize, m: &mut usize, samp: &Name) {
    let i = if *m == 0 { 0 } else { binfnd(sa, *m, samp) };
    if i != 0 {
        return;
    }
    *m += 1;
    if *m > mm {
        return;
    }
    sa[*m] = *samp;
    sa[1..=*m].sort_unstable();
}

// ---------------------------------------------------------------------------
// Sorting.  The Fortran routines are deterministic heapsorts:
//   * sort10/sort30/isort/sort order plain keys ascending; equal keys are
//     interchangeable, so Rust's sort_unstable gives identical results.
//   * sort2 orders by (dict, type) lexicographically ascending - including
//     its tie-breaking rules - so an ascending sort on the pair is exact.
// Fortran real comparison treats -0.0 == 0.0; partial_cmp does the same.
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself. `sort30` sorts the big arrays of (site, species, time)
/// or (neighbour, target, weight) triples built by `make_rec30`, so that records sharing the same
/// first field(s) end up next to each other and can be processed group by group. `sort_real` sorts
/// a plain array of numbers — used, for example, to find the percentile of estimated recording
/// effort across all sites at the end of a `frescalo` run.
pub fn sort30(dict: &mut [Rec30], n: usize) {
    dict[1..=n].sort_unstable();
}

pub fn sort_real(dict: &mut [f32], n: usize) {
    dict[1..=n].sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
}

/// Plumbing, not part of the method itself. `isort` sorts a plain array of integers, used in
/// `neighsim` to bring together the sites that share a given species so their pairwise similarity
/// can be tallied. `sort2` sorts two parallel arrays together by the values in the first one —
/// this is how `sampdist` and `neighsim` rank a site's neighbours by distance or similarity while
/// keeping track of which neighbour each ranked value belongs to, and how `frescalo`'s `fresca`
/// ranks species by rescaled frequency.
pub fn isort(dict: &mut [i32], n: usize) {
    dict[1..=n].sort_unstable();
}

pub fn sort2(dict: &mut [f32], typ: &mut [i32], n: usize) {
    let mut pairs: Vec<(f32, i32)> = (1..=n).map(|i| (dict[i], typ[i])).collect();
    pairs.sort_unstable_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then_with(|| a.1.cmp(&b.1))
    });
    for (k, (d, t)) in pairs.into_iter().enumerate() {
        dict[1 + k] = d;
        typ[1 + k] = t;
    }
}
