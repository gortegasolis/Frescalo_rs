//! Golden-output regression tests.
//!
//! The fixture in `tests/fixtures/` is synthetic (see `generate.py`). The golden
//! files in `tests/fixtures/golden/` were produced by a native gfortran build of
//! the original Fortran programs, and the Rust port must reproduce them byte for
//! byte. Any change to operation order, formatting or sorting shows up here.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const INPUTS: &[&str] = &[
    "Samp_locations.txt",
    "Samp_lonlat.txt",
    "Training.txt",
    "Test.txt",
    "NotBench.txt",
];

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures")
}

/// Fresh working directory holding a copy of the fixture inputs.
fn workdir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("frescalo_golden_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    for f in INPUTS {
        fs::copy(fixtures().join(f), dir.join(f)).unwrap();
    }
    dir
}

/// Run one of the binaries in `dir` with the given arguments and stdin.
fn run(bin: &str, dir: &Path, args: &[&str], stdin: &str) -> std::process::Output {
    let exe = match bin {
        "sampdist" => env!("CARGO_BIN_EXE_sampdist"),
        "neighsim" => env!("CARGO_BIN_EXE_neighsim"),
        "frescalo" => env!("CARGO_BIN_EXE_frescalo"),
        _ => unreachable!(),
    };
    let mut child = Command::new(exe)
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(stdin.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "{} failed: {}",
        bin,
        String::from_utf8_lossy(&out.stderr)
    );
    out
}

fn assert_same(dir: &Path, file: &str, golden: &str) {
    let got = fs::read(dir.join(file)).unwrap_or_else(|_| panic!("{} not written", file));
    let want = fs::read(fixtures().join("golden").join(golden).join(file)).unwrap();
    if got != want {
        let g = String::from_utf8_lossy(&got);
        let w = String::from_utf8_lossy(&want);
        let line = g
            .lines()
            .zip(w.lines())
            .position(|(a, b)| a != b)
            .map(|i| i + 1)
            .unwrap_or(g.lines().count().min(w.lines().count()) + 1);
        panic!("{} differs from golden/{} (first difference at line {})", file, golden, line);
    }
}

/// Run the whole pipeline through the interactive prompts, as
/// `verify_against_fortran.sh` does.
fn prompted_pipeline(dir: &Path) {
    run("sampdist", dir, &[], "Samp_locations.txt\ndist.txt\n50\n\n");
    run("neighsim", dir, &[], "Training.txt\ndist.txt\nsim.txt\nweights.txt\n25\n\n");
    run(
        "frescalo",
        dir,
        &[],
        "log.txt\nTest.txt\nweights.txt\n\nsamples.txt\nfrequencies.txt\ntrends.txt\n\n\n\n",
    );
}

#[test]
fn prompted_pipeline_matches_fortran() {
    let dir = workdir("prompted");
    prompted_pipeline(&dir);
    for f in ["dist.txt", "sim.txt", "weights.txt", "log.txt", "samples.txt", "frequencies.txt", "trends.txt"] {
        assert_same(&dir, f, "default");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn prompted_frescalo_with_exclusions_and_parameters_matches_fortran() {
    let dir = workdir("prompted_nondefault");
    fs::copy(fixtures().join("golden/default/weights.txt"), dir.join("weights.txt")).unwrap();
    run(
        "frescalo",
        &dir,
        &[],
        "log.txt\nTest.txt\nweights.txt\nNotBench.txt\nsamples.txt\nfrequencies.txt\ntrends.txt\n0.80\n0.15\n\n",
    );
    for f in ["log.txt", "samples.txt", "frequencies.txt", "trends.txt"] {
        assert_same(&dir, f, "nondefault");
    }
    let _ = fs::remove_dir_all(&dir);
}

/// Same pipeline driven entirely by options: must equal the prompted run.
#[test]
fn option_mode_pipeline_matches_fortran() {
    let dir = workdir("options");
    run(
        "sampdist",
        &dir,
        &["--locations", "Samp_locations.txt", "--output", "dist.txt", "--neighbours", "50", "--no-hold"],
        "",
    );
    run(
        "neighsim",
        &dir,
        &[
            "--training=Training.txt",
            "--distances=dist.txt",
            "--similarity-out=sim.txt",
            "--weights-out=weights.txt",
            "--neighbours=25",
            "--no-hold",
        ],
        "",
    );
    run(
        "frescalo",
        &dir,
        &[
            "--log", "log.txt",
            "--occurrences", "Test.txt",
            "--weights", "weights.txt",
            "--samples-out", "samples.txt",
            "--frequencies-out", "frequencies.txt",
            "--trends-out", "trends.txt",
            "--no-hold",
        ],
        "",
    );
    for f in ["dist.txt", "sim.txt", "weights.txt", "log.txt", "samples.txt", "frequencies.txt", "trends.txt"] {
        assert_same(&dir, f, "default");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn option_mode_frescalo_with_exclusions_and_parameters_matches_fortran() {
    let dir = workdir("options_nondefault");
    fs::copy(fixtures().join("golden/default/weights.txt"), dir.join("weights.txt")).unwrap();
    run(
        "frescalo",
        &dir,
        &[
            "--log", "log.txt",
            "--occurrences", "Test.txt",
            "--weights", "weights.txt",
            "--exclusions", "NotBench.txt",
            "--samples-out", "samples.txt",
            "--frequencies-out", "frequencies.txt",
            "--trends-out", "trends.txt",
            "--phi", "0.80",
            "--benchmark-limit", "0.15",
            "--no-hold",
        ],
        "",
    );
    for f in ["log.txt", "samples.txt", "frequencies.txt", "trends.txt"] {
        assert_same(&dir, f, "nondefault");
    }
    let _ = fs::remove_dir_all(&dir);
}

/// Options for some settings, prompts for the missing file names.
#[test]
fn mixed_options_and_prompts_match_fortran() {
    let dir = workdir("mixed");
    // output file prompted; the final <RETURN> answers the exit pause
    run("sampdist", &dir, &["--locations", "Samp_locations.txt", "--neighbours", "50"], "dist.txt\n\n");
    assert_same(&dir, "dist.txt", "default");
    let _ = fs::remove_dir_all(&dir);
}

/// A blank answer at the neighbours prompt takes the default (200, 100).
#[test]
fn blank_neighbours_answer_takes_default() {
    let dir = workdir("blank_neighbours");
    run("sampdist", &dir, &[], "Samp_locations.txt\ndist.txt\n\n\n");
    let dist = fs::read_to_string(dir.join("dist.txt")).unwrap();
    // 100 sites, fewer than the default 200: every site lists all 100
    assert_eq!(dist.lines().count(), 100 * 100);
    run("neighsim", &dir, &[], "Training.txt\ndist.txt\nsim.txt\nweights.txt\n\n\n");
    let weights = fs::read_to_string(dir.join("weights.txt")).unwrap();
    for line in weights.lines() {
        // columns: 2a10, 3f7.4, then K (i6) and N (i6)
        assert_eq!(line[41..47].trim(), "100", "{}", line);
        assert_eq!(line[47..53].trim(), "100", "{}", line);
    }
    let _ = fs::remove_dir_all(&dir);
}

/// Parse a sampdist row into (site, neighbour, rank, distance).
fn dist_row(line: &str) -> (String, String, usize, String) {
    (
        line[0..10].trim().to_string(),
        line[10..20].trim().to_string(),
        line[20..25].trim().parse().unwrap(),
        line[25..].trim().to_string(),
    )
}

/// Geodesic distances on longitude/latitude, across the 180° meridian.
#[test]
fn geodesic_distances_on_lon_lat() {
    let dir = workdir("geodesic");
    run(
        "sampdist",
        &dir,
        &["--locations", "Samp_lonlat.txt", "--output", "dist.txt", "--neighbours", "10", "--distance", "geodesic", "--no-hold"],
        "",
    );
    let text = fs::read_to_string(dir.join("dist.txt")).unwrap();
    let rows: Vec<_> = text.lines().map(dist_row).collect();
    assert_eq!(rows.len(), 100 * 10);
    let mut prev = -1.0f64;
    for (site, neighbour, rank, d) in &rows {
        assert!(!d.contains('*'), "overflowed distance: {}", d);
        let d: f64 = d.trim_end_matches('.').parse().unwrap();
        if *rank == 1 {
            assert_eq!(site, neighbour);
            assert_eq!(d, 0.0);
        } else {
            assert!(d >= prev, "distances not increasing for {}", site);
        }
        prev = d;
    }
    // S005 (179.0 W) and S004 (179.5 E) are 1.5 degrees apart across the
    // antimeridian (~137 km at 35 N): S004 is one of S005's nearest.
    let s004 = rows
        .iter()
        .find(|r| r.0 == "S005" && r.1 == "S004")
        .expect("S004 must be among S005's 10 nearest sites");
    assert!(s004.2 <= 3, "S004 ranked {}", s004.2);
    let km: f64 = s004.3.trim_end_matches('.').parse().unwrap();
    assert!((km - 137.0).abs() < 2.0, "S005-S004 distance {} km", km);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn geodesic_rejects_out_of_range_coordinates() {
    let dir = workdir("geodesic_bad");
    let out = Command::new(env!("CARGO_BIN_EXE_sampdist"))
        .args(["--locations", "Samp_locations.txt", "--output", "dist.txt", "--distance", "geodesic", "--no-hold"])
        .current_dir(&dir)
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("outside"));
    let _ = fs::remove_dir_all(&dir);
}
