//! Stage 3 (`frescalo`): a time factor for every species and period.
//!
//! [`hill`] holds the estimator; [`write_trends`] runs it over the records
//! (ordered by species, period, site) and writes the trends file, padding
//! with zero rows so that every species has one row per period.

pub mod hill;

use crate::dataset::Dataset;
use crate::fortran::{binfnd, cout, rec30_field, Arr2, Rec};
use crate::output::{self, write_trend, write_zero_trend};
use hill::tfcalc;
use std::io::Write;

/// Estimate and write the time factors. `ffij` holds the standardised local
/// frequencies, `sampef` the sampling intensities and `lendat` the number of
/// records per species and period.
#[allow(unused_assignments)] // jx/iitx updates mirror the Fortran control flow
pub fn write_trends<W: Write>(ds: &Dataset, ffij: &Arr2<f32>, sampef: &Arr2<f32>, lendat: &Arr2<i32>, out: &mut W) {
    let (m, n, nt) = (ds.m, ds.n, ds.nt);
    let (sp, tim, sa) = (&ds.sp, &ds.tim, &ds.sa);
    let mut iocc = vec![0i32; m + 2];
    let mut smpint = vec![0.0f32; m + 2];
    let mut fff = vec![0.0f32; m + 2];

    let mut idtji = 0usize;
    // Markers to ensure that time periods at which species were not recorded
    // are also written.
    let mut jx = 1usize;
    let mut iitx = 1usize;

    output::writeln(out, output::TRENDS_HEADER);

    // label 130 loop
    'main: loop {
        for i in 1..=m {
            iocc[i] = 0;
        }
        idtji += 1;
        if idtji > ds.nrec {
            break 'main;
        }
        let dtji = ds.records[idtji];
        let j = binfnd(sp, n, &rec30_field(&dtji, 0));
        let iit = binfnd(tim, nt, &rec30_field(&dtji, 1));
        for i in 1..=m {
            smpint[i] = sampef.at(i, iit);
            fff[i] = ffij.at(i, j);
        }
        let i = binfnd(sa, m, &rec30_field(&dtji, 2));
        iocc[i] = 1;
        for _ in 1..lendat.at(j, iit) {
            idtji += 1;
            if idtji > ds.nrec {
                break 'main;
            }
            let i = binfnd(sa, m, &rec30_field(&ds.records[idtji], 2));
            iocc[i] = 1;
        }
        let r = tfcalc(&iocc, &smpint, &fff, m);
        // first pad out the output with times when the species was not recorded
        if sp[jx] < sp[j] {
            for iiit in iitx..=nt {
                write_zero_trend(out, &sp[jx], &tim[iiit]);
            }
            for jj in (jx + 1)..j {
                for iiit in 1..=nt {
                    write_zero_trend(out, &sp[jj], &tim[iiit]);
                }
            }
            jx = j;
            iitx = 1;
        }
        if tim[iitx] < tim[iit] {
            for iiit in iitx..iit {
                write_zero_trend(out, &sp[j], &tim[iiit]);
            }
            iitx = iit;
        }
        write_trend(out, &sp[j], &tim[iit], r.tf, r.sd, r.jtot, r.sptot, r.esttot, r.ic1, r.ic2);
        if j % 10 == 0 {
            let mut rec = Rec::new();
            rec.name(&sp[j]).x(1).name(&tim[iit]).f(r.tf, 8, 3).f(r.sd, 7, 3).i(j as i64, 7);
            cout(&String::from_utf8_lossy(&rec.buf));
        }
        jx = j;
        iitx = iit + 1;
    }

    // label 140
    for iiit in iitx..=nt {
        write_zero_trend(out, &sp[jx], &tim[iiit]);
    }
}
