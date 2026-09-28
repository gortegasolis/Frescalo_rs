//! Rust port of Frescalo_1.f:
//! FRESCALO - Trend_analysis using local frequencies
//! written by Mark Hill, January-June 2011
//!
//! Calculates sampling effort multipliers and rescaled species probabilities
//! for each location ([`frescalo::effort`]), and time factors for each
//! species ([`frescalo::trend`]), following Hill (2012) Methods in Ecology
//! and Evolution 3: 195-205.

use frescalo::cli::{fatal, Args, OptSpec, Session, NO_HOLD};
use frescalo::config::{limits, FrescaloParams, BENCHMARK_LIMIT_RANGE, DEFAULT_BENCHMARK_LIMIT, DEFAULT_PHI, PHI_RANGE};
use frescalo::dataset::{read_exclusions, Dataset};
use frescalo::effort;
use frescalo::fortran::{blank_name, cout, ffmt, ifmt, read_f8_4, DataReader, LogWriter, Rec};
use frescalo::output::writeln;
use frescalo::trend::write_trends;
use std::fs::File;
use std::io::{BufWriter, Write};

const PROGRAM: &str = "frescalo";

const OPTIONS: &[OptSpec] = &[
    OptSpec {
        name: "log",
        value: Some("FILE"),
        help: "Log file (must not exist)",
    },
    OptSpec {
        name: "occurrences",
        value: Some("FILE"),
        help: "Occurrence records [site species period]",
    },
    OptSpec {
        name: "weights",
        value: Some("FILE"),
        help: "Neighbourhood weights from neighsim",
    },
    OptSpec {
        name: "exclusions",
        value: Some("FILE"),
        help: "Species to exclude from benchmarks, one per line [default: none]",
    },
    OptSpec {
        name: "samples-out",
        value: Some("FILE"),
        help: "Sample statistics output file (must not exist)",
    },
    OptSpec {
        name: "frequencies-out",
        value: Some("FILE"),
        help: "Rescaled frequency output file (must not exist)",
    },
    OptSpec {
        name: "trends-out",
        value: Some("FILE"),
        help: "Trend output file (must not exist)",
    },
    OptSpec {
        name: "phi",
        value: Some("PHI"),
        help: "Target local frequency, 0.50 to 0.95 [default: 0.74]",
    },
    OptSpec {
        name: "benchmark-limit",
        value: Some("R"),
        help: "Benchmark limit, 0.08 to 0.5 [default: 0.2703]",
    },
    NO_HOLD,
];

/// Write a line to the console and to the log.
fn both(log: &mut LogWriter, s: &str) {
    cout(s);
    writeln(log, s);
}

/// Report a fatal data error on the console and in the log, then stop.
fn fail(session: &mut Session, log: &mut LogWriter, msg: &str) -> ! {
    both(log, msg);
    log.flush().unwrap();
    session.fail();
}

/// Open an input file (option or prompt), echoing prompt and name to the log.
fn input(session: &mut Session, log: &mut LogWriter, given: Option<&str>, prompt: &str) -> File {
    writeln(log, prompt);
    let (name, f) = session.input_file(given, &[prompt]);
    writeln(log, &format!("{:<20}", name));
    f
}

/// Create an output file (option or prompt), echoing prompt and name to the log.
fn output(session: &mut Session, log: &mut LogWriter, given: Option<&str>, prompt: &str) -> BufWriter<File> {
    writeln(log, prompt);
    let (name, f) = session.output_file(given, &[prompt]);
    writeln(log, &format!("{:<20}", name));
    BufWriter::new(f)
}

/// Φ or R*: from the option (checked against `range`), the default in option
/// mode, or the original prompt loop (blank or 0 = default). The prompt line
/// always goes to the log, so option-mode logs match prompted ones.
fn parameter(
    session: &mut Session,
    log: &mut LogWriter,
    args: &Args,
    option: &str,
    prompt: &str,
    default: f32,
    range: (f32, f32),
    range_error: &str,
) -> f32 {
    if args.option_mode() {
        writeln(log, prompt);
        let v = args.parse_value::<f32>(PROGRAM, option).unwrap_or(default);
        let v = if v == 0.0 { default } else { v };
        if v > range.1 || v < range.0 {
            fatal(PROGRAM, &format!("--{} {}", option, range_error.trim_start_matches(" ***ERROR*** ")));
        }
        return v;
    }
    let mut value = default;
    loop {
        both(log, prompt);
        if let Some(v) = read_f8_4(session.stdin()) {
            value = v;
        }
        if value == 0.0 {
            value = default;
        }
        if value > range.1 || value < range.0 {
            cout(range_error);
            continue;
        }
        return value;
    }
}

fn main() {
    let args = Args::parse(PROGRAM, "FRESCALO - Trend_analysis using local frequencies", OPTIONS);
    let mut session = Session::new(PROGRAM, args.has("no-hold"));

    // Banner (format 1998; the leading // produces two blank records)
    let limits_lines = [
        " Input limits:".to_string(),
        format!("    Number of samples = {}", ifmt(limits::FRESCALO_SITES as i64, 7)),
        format!("    Number of species = {}", ifmt(limits::FRESCALO_SPECIES as i64, 7)),
        format!("    Number of time periods = {}", ifmt(limits::FRESCALO_PERIODS as i64, 7)),
        format!("    Number of observations = {}", ifmt(limits::FRESCALO_RECORDS as i64, 7)),
        format!("    Number of neighbourhood weights = {}", ifmt(limits::FRESCALO_WEIGHTS as i64, 7)),
        String::new(),
    ];
    cout("");
    cout("");
    cout(" FRESCALO - Trend_analysis using local frequencies");
    cout(" written by Mark Hill, January-June 2011");
    cout("");
    for l in &limits_lines {
        cout(l);
    }

    let log_prompt = " Type name of log file ...";
    let (fileou, flog) = session.output_file(args.get("log"), &[log_prompt]);
    let mut unit10 = BufWriter::new(flog);

    // Log banner (format 1999)
    writeln(&mut unit10, " Log file for FRESCALO");
    writeln(&mut unit10, "");
    for l in &limits_lines {
        writeln(&mut unit10, l);
    }
    writeln(&mut unit10, log_prompt);
    writeln(&mut unit10, &format!("{:<20}", fileou)); // 1000 format(a20)

    let focc = input(
        &mut session,
        &mut unit10,
        args.get("occurrences"),
        " Type occurrence input file [sample species time] ...",
    );
    let mut occ_reader = DataReader::new(focc);
    let fwgt = input(
        &mut session,
        &mut unit10,
        args.get("weights"),
        " Type neighbourhood weight input file ...",
    );
    let mut wgt_reader = DataReader::new(fwgt);

    // Benchmark exclusions
    let excl_prompt = [
        " Type file with species to exclude from benchmarks ",
        "     or press <RETURN> if no exclusions...",
    ];
    for p in excl_prompt {
        writeln(&mut unit10, p);
    }
    let excl_file = match args.get("exclusions") {
        Some(name) => match File::open(name) {
            Ok(f) => Some((name.to_string(), f)),
            Err(_) => fatal(PROGRAM, &format!("exclusion file '{}' does not exist", name)),
        },
        None if args.option_mode() => None,
        None => {
            for p in excl_prompt {
                cout(p);
            }
            // label 20: on open error the original silently re-reads,
            // without reprinting the prompt or any error message.
            loop {
                let name = session.read_a20();
                if name.is_empty() {
                    break None;
                }
                if let Ok(f) = File::open(&name) {
                    break Some((name, f));
                }
            }
        }
    };
    let bnchx = match excl_file {
        None => {
            writeln(&mut unit10, "<No exclusions>");
            vec![blank_name()]
        }
        Some((name, f)) => {
            writeln(&mut unit10, &format!("{:<20}", name));
            read_exclusions(f, limits::FRESCALO_SPECIES)
        }
    };
    let nbnchx = bnchx.len() - 1;

    let mut unit7 = output(
        &mut session,
        &mut unit10,
        args.get("samples-out"),
        " Type name of sample stats output file...",
    );
    let mut unit8 = output(
        &mut session,
        &mut unit10,
        args.get("frequencies-out"),
        " Type name of rescaled frequency file...",
    );
    let mut unit9 = output(
        &mut session,
        &mut unit10,
        args.get("trends-out"),
        " Type name of trend output file...",
    );

    // Target value of phi (label 40)
    let phi = parameter(
        &mut session,
        &mut unit10,
        &args,
        "phi",
        &format!(
            " Type target value of local frequency phi (default={})...",
            ffmt(DEFAULT_PHI, 4, 2)
        ),
        DEFAULT_PHI,
        PHI_RANGE,
        " ***ERROR*** Outside range 0.50 to 0.95",
    );
    // 2011 format has a trailing slash: blank line follows
    both(&mut unit10, &format!(" Target value is {}", ffmt(phi, 4, 2)));
    both(&mut unit10, "");

    // Benchmark limit (label 46)
    let blim = parameter(
        &mut session,
        &mut unit10,
        &args,
        "benchmark-limit",
        &format!(" Type value of Benchmark Limit (default={})...", ffmt(DEFAULT_BENCHMARK_LIMIT, 4, 2)),
        DEFAULT_BENCHMARK_LIMIT,
        BENCHMARK_LIMIT_RANGE,
        " ***ERROR*** Outside range 0.08 to 0.5",
    );
    // 2014 format has a trailing slash: blank line follows
    both(&mut unit10, &format!(" Benchmark limit is {}", ffmt(blim, 4, 2)));
    both(&mut unit10, "");

    let params = FrescaloParams {
        phi,
        benchmark_limit: blim,
        ..FrescaloParams::default()
    };

    // ------------------------------------------------------------------
    // Read in data
    // ------------------------------------------------------------------
    let mut ds = Dataset::new();
    if let Err(msg) = ds.read_weights(&mut wgt_reader) {
        fail(&mut session, &mut unit10, &msg);
    }
    if let Err(msg) = ds.read_occurrences(&mut occ_reader, &mut unit10) {
        fail(&mut session, &mut unit10, &msg);
    }

    // 2505 / 2506: report actual numbers
    let summary = ds.summary_lines(nbnchx);
    for l in &summary {
        cout(l);
    }
    for l in &summary {
        writeln(&mut unit10, l);
    }
    if nbnchx > 0 {
        writeln(&mut unit10, " Benchmark exclusions");
        for name in bnchx.iter().skip(1) {
            let mut r = Rec::new();
            r.x(4).name(name);
            r.writeln(&mut unit10);
        }
    }

    ds.sort_records();

    // Local frequencies, then effort standardisation site by site
    let totals = effort::weight_totals(&ds);
    let (idat, iitot) = effort::records_per_site(&ds);
    let mut lf = effort::local_frequencies(&ds, &idat, &totals);
    let bwght = effort::benchmark_weights(&ds, &bnchx);
    let mut st = effort::standardise_sites(&ds, &mut lf, &iitot, &totals, &bwght, &params, &mut unit7, &mut unit8);

    // Sampling intensity per site and period, then time factors
    ds.reorder_by_species_time();
    let (sampef, lendat) = effort::sampling_intensity(&ds, &st, &bwght);
    write_trends(&ds, &lf.ffij, &sampef, &lendat, &mut unit9);

    unit7.flush().unwrap();
    unit8.flush().unwrap();
    unit9.flush().unwrap();

    // Finally test whether given value of phi appears to be unrealistically low
    // (2513: leading // produces two blank records)
    let phi985 = effort::phi_in_985_percentile(&mut st.phi_in, ds.m);
    for _ in 0..2 {
        writeln(&mut unit10, "");
        cout("");
    }
    writeln(&mut unit10, &format!(" 98.5 percentile of input phi {}", ffmt(phi985, 5, 2)));
    writeln(&mut unit10, &format!(" Target value of phi          {}", ffmt(phi, 5, 2)));
    cout(&format!(" 98.5 percentile of input phi {}", ffmt(phi985, 5, 2)));
    cout(&format!(" Target value of phi          {}", ffmt(phi, 5, 2)));
    if phi < phi985 {
        // 2514: /// then text, /, text, /
        for _ in 0..3 {
            both(&mut unit10, "");
        }
        both(&mut unit10, " *** BEWARE *** ");
        both(&mut unit10, "");
        both(&mut unit10, " Target value of phi may be too small");
        both(&mut unit10, "");
    }
    // 2503: leading and trailing slashes
    writeln(&mut unit10, "");
    writeln(&mut unit10, " Calculation reached completion");
    writeln(&mut unit10, "");
    cout("");
    cout(" Calculation reached completion");
    cout("");
    unit10.flush().unwrap();

    session.finish();
}
