//! Writers for every output file, one function per line layout.
//!
//! Each layout is the original Fortran `format` statement (quoted in the
//! comments), so column widths, rounding and overflow asterisks are exactly
//! those of the originals.

use crate::fortran::{Name, Rec};
use std::io::Write;

/// Write `s` followed by a newline (log messages and table headers).
pub fn writeln<W: Write>(w: &mut W, s: &str) {
    w.write_all(s.as_bytes()).unwrap();
    w.write_all(b"\n").unwrap();
}

// ---------------------------------------------------------------------------
// sampdist
// ---------------------------------------------------------------------------

/// One `sampdist` row: site, neighbour, rank, distance.
/// `2030 format(2a10,i5,1x,f6.0)`
pub fn write_distance_row<W: Write>(w: &mut W, site: &Name, neighbour: &Name, rank: i64, dist: f32) {
    let mut r = Rec::new();
    r.name(site).name(neighbour).i(rank, 5).x(1).f(dist, 6, 0);
    r.writeln(w);
}

// ---------------------------------------------------------------------------
// neighsim
// ---------------------------------------------------------------------------

/// One similarity-file row: site, neighbour, similarity rank, Sørensen index.
/// `2030 format(2a10,i5,f6.3)`
pub fn write_similarity_row<W: Write>(w: &mut W, site: &Name, neighbour: &Name, rank: i64, sim: f32) {
    let mut r = Rec::new();
    r.name(site).name(neighbour).i(rank, 5).f(sim, 6, 3);
    r.writeln(w);
}

/// One weights-file row: site, neighbour, combined weight, similarity-rank
/// weight, distance-rank weight, K and N.
/// `2031 format(2a10,3f7.4,2i6)`
#[allow(clippy::too_many_arguments)]
pub fn write_weight_row<W: Write>(
    w: &mut W,
    site: &Name,
    neighbour: &Name,
    amult: f32,
    amult1: f32,
    amult2: f32,
    neigh: i32,
    neigh1: i32,
) {
    let mut r = Rec::new();
    r.name(site)
        .name(neighbour)
        .f(amult, 7, 4)
        .f(amult1, 7, 4)
        .f(amult2, 7, 4)
        .i(neigh as i64, 6)
        .i(neigh1 as i64, 6);
    r.writeln(w);
}

// ---------------------------------------------------------------------------
// frescalo
// ---------------------------------------------------------------------------

/// Header of the samples file (`2001`).
pub const SAMPLES_HEADER: &str = "Location  Loc_no  No_spp Phi_in  Alpha  Wgt_n2 Phi_out  Spnum_in Spnum_out Iter";

/// Header of the frequencies file (`2003`).
pub const FREQUENCIES_HEADER: &str = "Location   Species    Pres  Freq__  Freq_1 SD_Frq1  Rank  Rank_1";

/// Header of the trends file (`2060`).
pub const TRENDS_HEADER: &str = "Species__  Time______ TFactor St_Dev _Count ___spt ___est N>0.00 N>0.98";

/// One samples-file row (per site).
/// `2002 format(a10,2i7,f7.3,1x,f6.2,1x,f7.2,1x,f7.3,2f10.1,i5)`
#[allow(clippy::too_many_arguments)]
pub fn write_sample_row<W: Write>(
    w: &mut W,
    site: &Name,
    loc_no: usize,
    no_spp: i32,
    phi_in: f32,
    alpha: f32,
    wgt_n2: f32,
    phi_out: f32,
    spnum_in: f32,
    spnum_out: f32,
    iter: i32,
) {
    let mut r = Rec::new();
    r.name(site)
        .i(loc_no as i64, 7)
        .i(no_spp as i64, 7)
        .f(phi_in, 7, 3)
        .x(1)
        .f(alpha, 6, 2)
        .x(1)
        .f(wgt_n2, 7, 2)
        .x(1)
        .f(phi_out, 7, 3)
        .f(spnum_in, 10, 1)
        .f(spnum_out, 10, 1)
        .i(iter as i64, 5);
    r.writeln(w);
}

/// One frequencies-file row (per site and species).
/// `2004 format(a10,1x,a10,1x,i4,3f8.4,1x,i5,f8.3)`
#[allow(clippy::too_many_arguments)]
pub fn write_frequency_row<W: Write>(
    w: &mut W,
    site: &Name,
    species: &Name,
    pres: i32,
    freq: f32,
    freq1: f32,
    sd_frq1: f32,
    rank: usize,
    rank1: f32,
) {
    let mut r = Rec::new();
    r.name(site)
        .x(1)
        .name(species)
        .x(1)
        .i(pres as i64, 4)
        .f(freq, 8, 4)
        .f(freq1, 8, 4)
        .f(sd_frq1, 8, 4)
        .x(1)
        .i(rank as i64, 5)
        .f(rank1, 8, 3);
    r.writeln(w);
}

/// One trends-file row (per species and period).
/// `2050 format(a10,1x,a10,f8.3,f7.3,i7,2f7.1,2i7)`
#[allow(clippy::too_many_arguments)]
pub fn write_trend<W: Write>(
    w: &mut W,
    sp: &Name,
    tim: &Name,
    tf: f32,
    sd: f32,
    jtot: i32,
    spt: f32,
    est: f32,
    ic1: i32,
    ic2: i32,
) {
    let mut r = Rec::new();
    r.name(sp)
        .x(1)
        .name(tim)
        .f(tf, 8, 3)
        .f(sd, 7, 3)
        .i(jtot as i64, 7)
        .f(spt, 7, 1)
        .f(est, 7, 1)
        .i(ic1 as i64, 7)
        .i(ic2 as i64, 7);
    r.writeln(w);
}

/// A zero trends row: the species was not recorded in that period, but every
/// species still gets one row per period.
pub fn write_zero_trend<W: Write>(w: &mut W, sp: &Name, tim: &Name) {
    write_trend(w, sp, tim, 0.0, 0.0, 0, 0.0, 0.0, 0, 0);
}
