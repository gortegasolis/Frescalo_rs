//! Reading the `frescalo` inputs: neighbourhood weights, occurrence records
//! and the optional list of species excluded from the benchmarks.
//!
//! Sites are defined by the weights file: an occurrence record at a site that
//! has no neighbourhood is ignored (and noted in the log). Species and time
//! periods are whatever the occurrence file contains.

use crate::config::limits;
use crate::fortran::{
    addwrd, binfnd, blank_name, cout, ifmt, ld_i, ld_line, make_rec30, name_to_string, rec30_field, sort30, DataReader,
    Name, Rec, Rec30,
};
use std::fs::File;
use std::io::{self, BufRead, Write};

/// Everything `frescalo` reads. Lists are sorted and 1-based (element 0
/// unused), as in the original.
pub struct Dataset {
    /// Sites (from the weights file).
    pub sa: Vec<Name>,
    pub m: usize,
    /// Species (from the occurrence file).
    pub sp: Vec<Name>,
    pub n: usize,
    /// Time periods (from the occurrence file).
    pub tim: Vec<Name>,
    pub nt: usize,
    /// Occurrence records, `(site, species, period)` until
    /// [`Dataset::reorder_by_species_time`] turns them into
    /// `(species, period, site)`.
    pub records: Vec<Rec30>,
    pub nrec: usize,
    /// Weights stored as `(neighbour, target, weight)` and sorted, so that all
    /// the targets a site contributes to are together.
    pub weights: Vec<Rec30>,
    pub nwgt: usize,
    /// The weights reader's neighbour field as left after the last read; the
    /// original compares the first ignored occurrence site against it.
    weights_samp1: Name,
}

impl Dataset {
    pub fn new() -> Self {
        Dataset {
            sa: vec![blank_name(); limits::FRESCALO_SITES + 2],
            m: 0,
            sp: vec![blank_name(); limits::FRESCALO_SPECIES + 2],
            n: 0,
            tim: vec![blank_name(); limits::FRESCALO_PERIODS + 2],
            nt: 0,
            records: vec![[b' '; 30]; limits::FRESCALO_RECORDS + 2],
            nrec: 0,
            weights: vec![[b' '; 30]; limits::FRESCALO_WEIGHTS + 2],
            nwgt: 0,
            weights_samp1: blank_name(),
        }
    }

    /// Read the weights file (`target neighbour weight`), collecting the site
    /// list, then sort the weights. `Err` carries the limit message.
    pub fn read_weights(&mut self, reader: &mut DataReader) -> Result<(), String> {
        let mut samp = blank_name();
        let mut samp1 = blank_name();
        let mut weight = blank_name();
        loop {
            if self.nwgt == 0 {
                ld_line("Reading in smoothing weights from samples");
            }
            if !reader.getd(&mut samp, &mut samp1, &mut weight) {
                break;
            }
            self.nwgt += 1;
            if self.nwgt % 20000 == 0 {
                ld_line(&format!(
                    "Weights {}{}{}",
                    name_to_string(&samp),
                    name_to_string(&samp1),
                    ld_i(self.nwgt as i64)
                ));
            }
            if self.nwgt > limits::FRESCALO_WEIGHTS {
                return Err(format!(
                    " No of neighbourhood weights > maximum which is{}",
                    ifmt(limits::FRESCALO_WEIGHTS as i64, 5)
                ));
            }
            addwrd(&mut self.sa, &mut self.m, &samp);
            addwrd(&mut self.sa, &mut self.m, &samp1);
            self.weights[self.nwgt] = make_rec30(&samp1, &samp, &weight);
        }
        self.weights_samp1 = samp1;
        cout(" Sorting local frequency weights ...");
        sort30(&mut self.weights, self.nwgt);
        Ok(())
    }

    /// Read the occurrence file (`site species period`). Records at sites
    /// without a neighbourhood are skipped, with a note in `log` for each run
    /// of such records. `Err` carries the limit message.
    pub fn read_occurrences<L: Write>(&mut self, reader: &mut DataReader, log: &mut L) -> Result<(), String> {
        let mut samp = blank_name();
        // The original reuses the weights reader's neighbour field here as
        // "previous ignored site".
        let mut samp1 = self.weights_samp1;
        let mut spec = blank_name();
        let mut time = blank_name();
        while reader.getd(&mut samp, &mut spec, &mut time) {
            let i = binfnd(&self.sa, self.m, &samp);
            if i == 0 {
                if samp != samp1 {
                    // 2996
                    let mut r = Rec::new();
                    r.s("*** ")
                        .name(&samp)
                        .s(" location ignored - in species data but not listed in neighbourhood weights");
                    r.writeln(log);
                }
                samp1 = samp;
                continue;
            }
            addwrd(&mut self.sp, &mut self.n, &spec);
            addwrd(&mut self.tim, &mut self.nt, &time);
            if self.nt > limits::FRESCALO_PERIODS {
                return Err(format!(
                    " Number of time periods exceeds maximum which is{}",
                    ifmt(limits::FRESCALO_PERIODS as i64, 5)
                ));
            }
            self.nrec += 1;
            if self.nrec % 20000 == 0 {
                ld_line(&format!(
                    "{}{}{}{}",
                    name_to_string(&samp),
                    name_to_string(&spec),
                    name_to_string(&time),
                    ld_i(self.nrec as i64)
                ));
            }
            self.records[self.nrec] = make_rec30(&samp, &spec, &time);
        }
        Ok(())
    }

    /// Sort the occurrence records by (site, species, period).
    pub fn sort_records(&mut self) {
        cout(" Sorting main data ...");
        sort30(&mut self.records, self.nrec);
    }

    /// Reorder the occurrence records to (species, period, site) and sort
    /// them, ready for the time-factor stage.
    pub fn reorder_by_species_time(&mut self) {
        for idtji in 1..=self.nrec {
            let dtji = self.records[idtji];
            let f_samp = rec30_field(&dtji, 0);
            let f_spec = rec30_field(&dtji, 1);
            let f_time = rec30_field(&dtji, 2);
            if idtji % 20000 == 0 {
                ld_line(&format!(
                    "Re-ordering {}{}{}{}",
                    name_to_string(&f_samp),
                    name_to_string(&f_spec),
                    name_to_string(&f_time),
                    ld_i(idtji as i64)
                ));
            }
            self.records[idtji] = make_rec30(&f_spec, &f_time, &f_samp);
        }
        cout(" Now doing second sort of reordered main data ...");
        sort30(&mut self.records, self.nrec);
    }

    /// The "Actual numbers in data" summary lines (format 2505/2506).
    pub fn summary_lines(&self, nbnchx: usize) -> Vec<String> {
        vec![
            String::new(),
            " Actual numbers in data".to_string(),
            format!("    Number of samples      {}", ifmt(self.m as i64, 8)),
            format!("    Number of species      {}", ifmt(self.n as i64, 8)),
            format!("    Number of time periods {}", ifmt(self.nt as i64, 8)),
            format!("    Number of observations {}", ifmt(self.nrec as i64, 8)),
            format!("    Neighbourhood weights  {}", ifmt(self.nwgt as i64, 8)),
            format!("    Benchmark exclusions   {}", ifmt(nbnchx as i64, 8)),
            String::new(),
        ]
    }
}

impl Default for Dataset {
    fn default() -> Self {
        Self::new()
    }
}

/// Read the benchmark-exclusion file: one species per line, first 10
/// characters (Fortran `a10`), at most `max` species. Returns the list
/// 1-based (element 0 unused).
pub fn read_exclusions(f: File, max: usize) -> Vec<Name> {
    let mut out = vec![blank_name()];
    let mut lines = io::BufReader::new(f).lines();
    for _ in 1..=max {
        match lines.next() {
            Some(Ok(line)) => {
                let bytes = line.as_bytes();
                let mut len = bytes.len();
                if len > 0 && bytes[len - 1] == b'\r' {
                    len -= 1;
                }
                let mut nm = blank_name();
                let c = len.min(10);
                nm[..c].copy_from_slice(&bytes[..c]);
                out.push(nm);
            }
            _ => break,
        }
    }
    out
}
