//! Fortran-style formatted output (`Iw`, `Fw.d`, fixed-column records).

use super::text::{Name, BLANK};
use std::io::Write;

// ---------------------------------------------------------------------------
// Fortran-style formatted output.
// ---------------------------------------------------------------------------

/// Plumbing, not part of the method itself. Reproduces Fortran's classic fixed-width number
/// formatting (`Iw` for integers, `Fw.d` for decimals), right-aligned in a field of `w` characters
/// and replaced with `*`s if the number doesn't fit — the same formatting rule the original `.exe`
/// programs used, which is why the output files this port produces (e.g. `samples.txt`,
/// `trends.txt`) line up column-for-column identically to the originals, byte for byte (see the
/// project README's fidelity notes).
/// Integer edit descriptor Iw: right-justified, asterisks on overflow.
pub fn ifmt(v: i64, w: usize) -> String {
    let s = v.to_string();
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{:>width$}", s, width = w)
    }
}

/// Real edit descriptor Fw.d.
/// * rounds to d decimals (correct rounding of the exact binary value);
/// * Fw.0 keeps the trailing decimal point, as Fortran does ("  123.");
/// * if the field overflows, a leading zero is dropped ("0.74" -> ".74"),
///   and asterisks are printed if it still does not fit;
/// * Infinity/NaN follow gfortran conventions.
pub fn ffmt(v: f32, w: usize, d: usize) -> String {
    if v.is_nan() {
        return if w >= 3 {
            format!("{:>width$}", "NaN", width = w)
        } else {
            "*".repeat(w)
        };
    }
    if v.is_infinite() {
        let (full, short) = if v > 0.0 {
            ("Infinity", "Inf")
        } else {
            ("-Infinity", "-Inf")
        };
        if w >= full.len() {
            return format!("{:>width$}", full, width = w);
        }
        if w >= short.len() {
            return format!("{:>width$}", short, width = w);
        }
        return "*".repeat(w);
    }
    let mut s = format!("{:.*}", d, v);
    if d == 0 {
        s.push('.');
    }
    if s.len() > w {
        if let Some(rest) = s.strip_prefix("0.") {
            s = format!(".{}", rest);
        } else if let Some(rest) = s.strip_prefix("-0.") {
            s = format!("-.{}", rest);
        }
    }
    if s.len() > w {
        return "*".repeat(w);
    }
    format!("{:>width$}", s, width = w)
}

/// Record builder for formatted output lines (byte-oriented so that a10
/// name fields are reproduced exactly).
pub struct Rec {
    pub buf: Vec<u8>,
}

impl Rec {
    /// Plumbing, not part of the method itself. `Rec` builds up one output line field by field (a site
    /// name, then a number, then another number, ...), mimicking Fortran's fixed-column `WRITE`
    /// statements. `new` starts an empty line; `raw` appends already-formatted bytes verbatim (used
    /// when copying a 30-character record straight into a log message). See the `s`/`name`/`x`/`i`/`f`
    /// methods below for the accessors that actually add formatted fields.
    pub fn new() -> Self {
        Rec { buf: Vec::new() }
    }
    pub fn raw(&mut self, s: &[u8]) -> &mut Self {
        self.buf.extend_from_slice(s);
        self
    }
    /// Plumbing, not part of the method itself. These append one more field to a `Rec` output line
    /// under construction: `s` for raw text, `name` for a site/species/time code, `x` for `n` blank
    /// characters (column padding/spacing), `i` for a right-aligned whole number, `f` for a
    /// right-aligned decimal number. `writeln` finishes the line and writes it out. Chaining these
    /// calls is what produces each row of, for example, the sample-statistics or rescaled-frequency
    /// output tables — the Rust equivalent of the original's fixed-format `WRITE` statements.
    pub fn s(&mut self, s: &str) -> &mut Self {
        self.buf.extend_from_slice(s.as_bytes());
        self
    }
    pub fn name(&mut self, n: &Name) -> &mut Self {
        self.buf.extend_from_slice(n);
        self
    }
    pub fn x(&mut self, n: usize) -> &mut Self {
        self.buf.extend(std::iter::repeat(BLANK).take(n));
        self
    }
    pub fn i(&mut self, v: i64, w: usize) -> &mut Self {
        self.s(&ifmt(v, w));
        self
    }
    pub fn f(&mut self, v: f32, w: usize, d: usize) -> &mut Self {
        self.s(&ffmt(v, w, d));
        self
    }
    pub fn writeln<W: Write>(&self, w: &mut W) {
        w.write_all(&self.buf).unwrap();
        w.write_all(b"\n").unwrap();
    }
}
