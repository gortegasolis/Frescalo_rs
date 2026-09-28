//! Stage 1a (`sampdist`): each site's geographic neighbours, ranked closest
//! first.
//!
//! Given one line per site (`site x y`), every pairwise distance is computed
//! and, for each site, its *N* nearest sites are written out in order of
//! increasing distance, the site itself first. This is the spatial shortlist
//! from which `neighsim` later picks the floristically most similar sites
//! (Hill 2012 used N = 200 hectads).

use crate::distance::{Coord, DistanceEngine};
use crate::fortran::{addwrd, binfnd, blank_name, getnum, getnum_f64, ld_f, ld_i, ld_line, name_to_string, DataReader, Name};

/// Site codes (sorted, 1-based) and their coordinates.
pub struct Locations {
    pub names: Vec<Name>,
    pub coords: Vec<Coord>,
    pub m: usize,
}

/// Read the locations file: a first pass collects and sorts the site codes,
/// a second pass reads the coordinates (as the original does).
pub fn read_locations(reader: &mut DataReader, max_sites: usize) -> Locations {
    let mut names = vec![blank_name(); max_sites + 2];
    let mut coords = vec![Coord::default(); max_sites + 2];
    let mut m: usize = 0;

    let mut samp = blank_name();
    let mut east = blank_name();
    let mut north = blank_name();
    while reader.getd(&mut samp, &mut east, &mut north) {
        addwrd(&mut names, &mut m, &samp);
        if m % 100 == 0 {
            ld_line(&format!("{}  Sample  {}", ld_i(m as i64), name_to_string(&samp)));
        }
    }

    reader.rewind();

    while reader.getd(&mut samp, &mut east, &mut north) {
        let i = binfnd(&names, m, &samp);
        coords[i] = Coord {
            x32: getnum(&east),
            y32: getnum(&north),
            x64: getnum_f64(&east),
            y64: getnum_f64(&north),
        };
        if i % 100 == 0 {
            ld_line(&format!(
                "{}   {}{}{}",
                ld_i(i as i64),
                name_to_string(&names[i]),
                ld_f(coords[i].x32),
                ld_f(coords[i].y32)
            ));
        }
    }

    Locations { names, coords, m }
}

/// True if every coordinate pair lies within longitude/latitude bounds, a hint
/// that planar mode may have been chosen for data in degrees.
pub fn looks_like_degrees(locs: &Locations) -> bool {
    locs.m > 0
        && (1..=locs.m).all(|i| {
            let c = &locs.coords[i];
            (-180.0..=180.0).contains(&c.x64) && (-90.0..=90.0).contains(&c.y64)
        })
}

/// All sites ranked by distance from site `i1` (ties broken by site number),
/// into `ranked` as `(distance, site number)`. Equivalent to the original's
/// `sort2` on single-precision distances in planar mode.
pub fn rank_by_distance(i1: usize, locs: &Locations, engine: &DistanceEngine, ranked: &mut Vec<(f64, i32)>) {
    ranked.clear();
    for i2 in 1..=locs.m {
        ranked.push((engine.distance(&locs.coords[i1], &locs.coords[i2]), i2 as i32));
    }
    ranked.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then_with(|| a.1.cmp(&b.1)));
}
