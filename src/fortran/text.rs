//! Fixed-width text fields and the input-file reader (Fortran `getd`/`getnum`).

use std::fs::File;
use std::io::{BufRead, BufReader};

pub type Name = [u8; 10];
pub type Rec30 = [u8; 30];

pub const BLANK: u8 = b' ';

/// Plumbing, not part of the method itself. Site codes, species codes and time-period labels are
/// all handled as fixed-width, 10-character, blank-padded text fields throughout the pipeline
/// (matching the Fortran `character*10` fields of the original file formats). `blank_name` creates
/// an empty one of these fields; `name_to_string` converts one back to an ordinary trimmed string
/// for printing to the screen or a log file.
pub fn blank_name() -> Name {
    [BLANK; 10]
}

pub fn name_to_string(n: &Name) -> String {
    String::from_utf8_lossy(n).trim_end().to_string()
}

/// Plumbing, not part of the method itself. An occurrence record is really a triple — e.g. (site,
/// species, time period), or (neighbour site, target site, weight) — and `frescalo`/`neighsim`
/// sort huge numbers of these triples so that matching records end up adjacent (all records for
/// one site together, sorted by species; all a target site's weighted neighbours together, etc.).
/// `make_rec30` glues three 10-character fields into one 30-character record so they sort and
/// compare as a single unit; `rec30_field` pulls one of the three fields back out again.
pub fn make_rec30(a: &Name, b: &Name, c: &Name) -> Rec30 {
    let mut r = [BLANK; 30];
    r[0..10].copy_from_slice(a);
    r[10..20].copy_from_slice(b);
    r[20..30].copy_from_slice(c);
    r
}

/// field index: 0, 1 or 2 (Fortran d(1), d(2), d(3))
pub fn rec30_field(r: &Rec30, i: usize) -> Name {
    let mut n = [BLANK; 10];
    n.copy_from_slice(&r[i * 10..i * 10 + 10]);
    n
}

// ---------------------------------------------------------------------------
// getd: read one record and split into the first three blank-delimited words.
// Faithful port of the Fortran subroutine getd, including these quirks:
//   * words 1 and 2 are truncated to 10 characters, word 3 to 9 characters;
//   * the scan for the next word starts two positions after the end of the
//     previous word (k2+2), matching the original exactly;
//   * a blank line, or a line with fewer words than expected, returns with
//     the output fields UNCHANGED (the caller then re-processes stale data);
//   * end of file returns `false` (Fortran iend = 1).
// ---------------------------------------------------------------------------
pub struct DataReader {
    reader: BufReader<File>,
}

impl DataReader {
    /// Plumbing, not part of the method itself. `DataReader` reads the plain-text input files line by
    /// line and splits each line into (up to) three whitespace-separated fields — this is how a line
    /// like `NB13 Amblystegium 1990` in an occurrence file, or `NB13 ND04 0.108` in a neighbourhood
    /// weights file, gets turned into the three fixed-width text fields the rest of the program works
    /// with. `rewind` lets a file be read through twice (`sampdist` reads its location file once to
    /// collect all site codes, then again to read the coordinates); `getd` reads one line and returns
    /// `false` once the file is exhausted.
    pub fn new(f: File) -> Self {
        DataReader {
            reader: BufReader::new(f),
        }
    }

    /// Fortran rewind: seeking a BufReader discards its buffer.
    pub fn rewind(&mut self) {
        use std::io::Seek;
        let _ = self.reader.seek(std::io::SeekFrom::Start(0));
    }

    /// Returns false at end of file (iend = 1).
    pub fn getd(&mut self, w1: &mut Name, w2: &mut Name, w3: &mut Name) -> bool {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => return false, // EOF -> iend = 1
            Ok(_) => {}
            Err(_) => return false,
        }
        if line.ends_with('\n') {
            line.pop();
        }
        // Fortran read with format a80: take at most 80 bytes, blank-padded.
        // Windows CRLF files: strip the trailing CR, as the records seen by
        // the original Windows executables never contained it.
        let bytes = line.as_bytes();
        let mut len = bytes.len();
        if len > 0 && bytes[len - 1] == b'\r' {
            len -= 1;
        }
        let mut b = [BLANK; 80];
        let n = len.min(80);
        b[..n].copy_from_slice(&bytes[..n]);

        // k1 = first non-blank
        let mut k = 0usize;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        if k >= 80 {
            return true; // blank record: return with fields unchanged
        }
        let k1 = k;
        // k2 = last character of word 1
        let mut k = k1;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        if k >= 80 {
            return true; // no terminating blank: return unchanged
        }
        let k2 = k - 1;
        // k3 = first character of word 2 (scan starts at k2+2)
        let mut k = k2 + 2;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        if k >= 80 {
            return true;
        }
        let k3 = k;
        // k4 = last character of word 2
        let mut k = k3;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        if k >= 80 {
            return true;
        }
        let k4 = k - 1;
        // k5 = first character of word 3 (scan starts at k4+2)
        let mut k = k4 + 2;
        while k < 80 && b[k] == BLANK {
            k += 1;
        }
        let k5 = if k >= 80 { 79 } else { k }; // Fortran: if(k.gt.79) k=80
        // k6 = last character of word 3
        let mut k = k5;
        while k < 80 && b[k] != BLANK {
            k += 1;
        }
        let mut k6 = if k >= 80 { 79 } else { k - 1 };
        if k5 == 79 {
            k6 = 79;
        }

        let mut nw1 = [BLANK; 10];
        for k in k1..=k2 {
            if k < k1 + 10 {
                nw1[k - k1] = b[k];
            }
        }
        let mut nw2 = [BLANK; 10];
        for k in k3..=k4 {
            if k < k3 + 10 {
                nw2[k - k3] = b[k];
            }
        }
        let mut nw3 = [BLANK; 10];
        for k in k5..=k6 {
            if k < k5 + 9 {
                nw3[k - k5] = b[k];
            }
        }
        *w1 = nw1;
        *w2 = nw2;
        *w3 = nw3;
        true
    }
}

// ---------------------------------------------------------------------------
// getnum: parse a real number from a character*10 word (Fortran getnum).
// If the word has no decimal point before its first blank, a point is appended
// at the first blank position (this defeats the implied-decimals rule of the
// F10.4 descriptor).  Blanks are then ignored (BN mode) and the rest parsed;
// any parse error yields 0 (Fortran err=200 branch).
// ---------------------------------------------------------------------------
/// Converts a 10-character text field into an actual floating-point number. Its main job is
/// turning the neighbourhood-weight column written by `neighsim` — the `w_ii'` weight for one
/// neighbour pair, e.g. 0.108 in the paper's worked example (Hill 2012, p.197) — from text back
/// into a number `frescalo` can multiply and sum with. It's plumbing: the numeric parsing itself
/// isn't part of the method, only the value it recovers is.
pub fn getnum(weight: &Name) -> f32 {
    match getnum_text(weight) {
        Some(s) => s.parse::<f32>().unwrap_or(0.0),
        None => 0.0,
    }
}

/// Double-precision variant of [`getnum`], with the same text handling; used
/// for longitude/latitude in geodesic mode.
pub fn getnum_f64(weight: &Name) -> f64 {
    match getnum_text(weight) {
        Some(s) => s.parse::<f64>().unwrap_or(0.0),
        None => 0.0,
    }
}

/// The text `getnum` hands to the number parser (None for an empty field).
fn getnum_text(weight: &Name) -> Option<String> {
    let mut w = *weight;
    let mut idot = false;
    let mut kblank: Option<usize> = None;
    for k in 0..10 {
        if w[k] == b'.' {
            idot = true;
        }
        if w[k] == BLANK {
            kblank = Some(k);
            break;
        }
    }
    if !idot {
        if let Some(k) = kblank {
            w[k] = b'.';
        }
        // (If all 10 characters are non-blank with no dot, the original
        //  writes past w(10); we simply leave the word unchanged.)
    }
    let s: Vec<u8> = w.iter().filter(|&&c| c != BLANK).copied().collect();
    if s.is_empty() {
        return None;
    }
    Some(String::from_utf8_lossy(&s).into_owned())
}
