//! Console prompts, progress messages and the `filin`/`filout`/`hold` helpers.

use super::text::BLANK;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};

// ---------------------------------------------------------------------------
// Console (stdout) output.  Prompts are flushed so that they appear before
// input is read.  List-directed output (write(*,*)) is approximated: leading
// blank, character items concatenated, integers in a field of width 12.
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself: prints one line to the screen, flushed immediately.
/// Used for every interactive prompt and progress message the programs print while running.
pub fn cout(s: &str) {
    let mut o = io::stdout();
    let _ = o.write_all(s.as_bytes());
    let _ = o.write_all(b"\n");
    let _ = o.flush();
}

/// Plumbing, not part of the method itself: formats numbers the way Fortran's free-form
/// ("list-directed") console output does, and prints the periodic progress lines you see scroll by
/// while a large data set is being read, sorted, or rescaled — datasets in the paper's own example
/// ran to over 2 million bryophyte records (Hill 2012, p.196), so this kind of "still working..."
/// feedback matters in practice.
pub fn ld_i(v: i64) -> String {
    format!("{:12}", v)
}

pub fn ld_f(v: f32) -> String {
    format!(" {:>12}", format!("{}", v))
}

/// Approximate a Fortran list-directed write: a leading blank followed by the
/// (already concatenated) items.
pub fn ld_line(body: &str) {
    cout(&format!(" {}", body));
}

// ---------------------------------------------------------------------------
// Console input.  All console reads in the originals are formatted or
// list-directed reads from unit *; hitting EOF without an end= branch is a
// runtime error in Fortran, so we exit with an error message.
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself. `read_line` reads one line typed by the user (or
/// piped in), exiting with an error message if input runs out unexpectedly — matching the original
/// Fortran's behaviour on end-of-file. `read_a20` reads a line and keeps only the first 20
/// characters, used for every "type a file name" prompt (input files, output files, the optional
/// benchmark-exclusion file).
pub fn read_line(stdin: &mut impl BufRead) -> String {
    let mut s = String::new();
    match stdin.read_line(&mut s) {
        Ok(0) => {
            eprintln!("\nFortran runtime error: End of file");
            std::process::exit(2);
        }
        Ok(_) => {
            if s.ends_with('\n') {
                s.pop();
                if s.ends_with('\r') {
                    s.pop();
                }
            }
            s
        }
        Err(e) => {
            eprintln!("\nFortran runtime error: {}", e);
            std::process::exit(2);
        }
    }
}

/// read(*,1000) filein  with format a20: first 20 characters of the record,
/// trailing blanks trimmed for OPEN.
pub fn read_a20(stdin: &mut impl BufRead) -> String {
    let line = read_line(stdin);
    let bytes = line.as_bytes();
    let n = bytes.len().min(20);
    String::from_utf8_lossy(&bytes[..n]).trim_end().to_string()
}

/// Plumbing, not part of the method itself, but these two read the two parameters that matter most
/// to the method's results: `read_f8_4` reads the target local frequency Φ and the benchmark limit
/// R* typed at the `frescalo` prompts (pressing return keeps the default, Φ=0.74, R*=0.2703 — see
/// `frescalo/main` for what those mean). `read_int_listdirected` reads a plain whole number — used
/// by `sampdist` and `neighsim` for "how many neighbours to include".
/// read(*,2010,err=...) x  with format f8.4 (BN blank mode):
/// blanks are ignored; an all-blank field reads as 0; without a decimal
/// point the field has 4 implied decimals.  Returns None on the err= branch.
pub fn read_f8_4(stdin: &mut impl BufRead) -> Option<f32> {
    let line = read_line(stdin);
    let bytes = line.as_bytes();
    let n = bytes.len().min(8);
    let field: Vec<u8> = bytes[..n].iter().filter(|&&c| c != BLANK).copied().collect();
    if field.is_empty() {
        return Some(0.0);
    }
    let s = String::from_utf8_lossy(&field).into_owned();
    if s.contains('.') {
        s.parse::<f32>().ok()
    } else {
        match s.parse::<i64>() {
            Ok(v) => Some(v as f32 / 10000.0),
            Err(_) => None,
        }
    }
}

/// read(*,*) neigh : list-directed integer read.
pub fn read_int_listdirected(stdin: &mut impl BufRead) -> i32 {
    let line = read_line(stdin);
    match line.split_whitespace().next() {
        Some(tok) => match tok.parse::<i32>() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("\nFortran runtime error: Bad integer for item 1 in list input");
                std::process::exit(2);
            }
        },
        None => {
            eprintln!("\nFortran runtime error: End of file");
            std::process::exit(2);
        }
    }
}

// ---------------------------------------------------------------------------
// filin / filout / hold
// ---------------------------------------------------------------------------
/// Plumbing, not part of the method itself. `filin`/`filout` implement the "type a file name"
/// prompts, re-prompting if the input file doesn't exist or the output file would overwrite an
/// existing one (the originals refuse to overwrite outputs, so re-running a program with the same
/// output name fails safely rather than silently clobbering a previous result). `hold` prints
/// "Press \<RETURN\> to exit" and waits — the modal pause the original console programs show once
/// a run has finished (or hit a fatal data-size limit) before the window closes.
pub fn filin(stdin: &mut impl BufRead) -> (String, File) {
    loop {
        let name = read_a20(stdin);
        match File::open(&name) {
            Ok(f) => return (name, f),
            Err(_) => {
                cout("");
                cout("  *** ERROR *** File does not exist");
                cout(" Type another name");
            }
        }
    }
}

pub fn filout(stdin: &mut impl BufRead) -> (String, File) {
    loop {
        let name = read_a20(stdin);
        match OpenOptions::new().write(true).create_new(true).open(&name) {
            Ok(f) => return (name, f),
            Err(_) => {
                cout("");
                cout("  *** ERROR *** File already exists");
                cout(" Type another name");
            }
        }
    }
}

/// Fortran hold: prompt and wait for `<RETURN>`, then stop.
pub fn hold(stdin: &mut impl BufRead) -> ! {
    let mut o = io::stdout();
    let _ = o.write_all(b"\n\nPress <RETURN> to exit\n\n\n----------------------\n");
    let _ = o.flush();
    let mut s = String::new();
    let _ = stdin.read_line(&mut s);
    std::process::exit(0);
}
