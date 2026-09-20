# Frescalo2Rust

A faithful, behaviour-preserving Rust port of Mark Hill's FRESCALO suite (Hill, 2011, *Methods in Ecology and Evolution* 2: 502–512; see `Original_frescalo/` for the original Fortran sources and the paper).

Three programs are ported, as three binaries:

| Binary | Original | Purpose |
|----------------------|--------------------------------|------------------|
| `sampdist` | `sampdist_1.f` | Euclidean distances between sample locations; writes nearest neighbours per location. |
| `neighsim` | `neighsim_1.f` | Floristic similarity + physical proximity; writes neighbourhood weights. |
| `frescalo` | `Frescalo_1.f` | Sampling-effort multipliers, rescaled species frequencies, and species time factors. |

## Build and run

``` sh
cargo build --release
./target/release/sampdist      # interactive, exactly like the originals
./target/release/neighsim
./target/release/frescalo
```

The programs are interactive and prompt for file names and parameters on stdin, reproducing the originals' prompts; they can be driven with piped input, e.g.:

``` sh
printf 'log.txt\nTest.txt\nweights.txt\n\nsamples.txt\nfrequencies.txt\ntrends.txt\n\n\n\n' \
  | ./target/release/frescalo
```

(Inputs: log file, occurrence file, weights file, benchmark-exclusion file \[blank = none\], three output files, target phi \[blank = 0.74\], benchmark limit \[blank = 0.2703\], and a final <RETURN>.)

## Documentation

Every function has a `///` doc comment (ported from the project's earlier Quarto book, which has been retired in favour of Rust's own tooling). Generate and open the HTML docs with:

``` sh
cargo doc --no-deps --open
```

## Development script

`./frescalo.sh` wraps the common build/test/release tasks; run `./frescalo.sh help` for the full list. Highlights:

``` sh
./frescalo.sh build     # cargo build --release
./frescalo.sh doc       # cargo doc --no-deps --open
./frescalo.sh verify    # run verify_against_fortran.sh (needs Original_frescalo/)
./frescalo.sh release   # tag a version and push it, triggering the release workflow
```

## Continuous integration

- `.github/workflows/ci.yml` builds and tests on every push/PR for Linux, Windows, and macOS, and uploads each platform's `sampdist`/`neighsim`/ `frescalo` binaries as workflow artifacts.
- `.github/workflows/release.yml` runs on `v*` tags: it builds release binaries for all three platforms, packages them into archives, and attaches them to a GitHub Release.
- `.github/workflows/docs.yml` builds `cargo doc` on pushes to `master` and publishes it to GitHub Pages.

### Compiling the binaries via GitHub Actions

No local Rust toolchain is needed to get `sampdist`, `neighsim`, and `frescalo` for Linux, Windows, or macOS — GitHub Actions builds all three platforms on every push and on tagged releases.

- **Workflow artifacts (any commit or PR)**: `ci.yml` runs automatically on every push to `main`/`master` and on pull requests, and can also be started manually from the repo's **Actions → CI → Run workflow** button (or `gh workflow run ci.yml`), since it declares `workflow_dispatch`. Open the finished run and download `frescalo-linux-x86_64`, `frescalo-windows-x86_64`, or `frescalo-macos-arm64` from the run's *Artifacts* section — each contains that platform's three binaries (`.exe` on Windows).

- **Tagged releases (packaged archives)**: pushing a `v*` tag triggers `release.yml`, which builds release binaries for all three platforms, packages each into a `frescalo-<tag>-<platform>.{tar.gz,zip}` archive (bundled with `README.md`), and attaches all three archives to a new GitHub Release. To cut one:

  ``` sh
  ./frescalo.sh release   # bumps/tags the version and pushes the tag
  # or manually:
  git tag v0.1.0
  git push origin v0.1.0
  ```

  Then grab the platform archive you need from the release's *Assets*.

Fortran verification (`verify_against_fortran.sh`) needs `Original_frescalo/`, which is gitignored and only available locally, so it is not part of CI — run it locally after checking out this repo alongside the original Fortran sources and data files.

## Verification

``` sh
./verify_against_fortran.sh
```

builds the original Fortran with gfortran and diffs every output file of the full pipeline (defaults, and a second run with `NotBench.txt` exclusions plus non-default phi/blim). **All output files are byte-identical** to the native Fortran build on this platform.

### Relationship to the Windows reference outputs

`Original_frescalo/` contains output files produced by the original Windows executables. The Rust port reproduces them exactly except for last-digit rounding in a small fraction of lines (`dist.txt` 0/80800, `samples.txt` 0/405, `trends.txt` 3/6521, `frequencies.txt` 583/185803, `sim.txt` 16/40400, `weights.txt` 12/35549). These differences are platform artifacts of the old Windows compiler: it evaluated intermediate expressions in x87 80-bit extended precision and rounded exact ties half-away-from-zero, whereas modern x86-64 (both gfortran and this port) uses strict IEEE single precision and round-half-to-even. A native gfortran build of the original sources produces output identical to this port, not to the Windows reference files.

## Fidelity notes

The port deliberately preserves the original numerics and quirks:

- **Single precision**: all Fortran `real` arithmetic is `f32`, with the same operation order, so results match a native Fortran build bit-for-bit.
- **Fixed-width text**: names are blank-padded 10-byte fields (9 bytes for the third word of a data line), compared byte-wise, and packed into 30-byte sort records, exactly like the Fortran `character` variables.
- **`getd` parser quirks**: the scan for the next word starts two characters after the end of the previous one; blank or short lines leave the fields unchanged (the caller then re-processes the previous values).
- **`getnum`**: a decimal point is appended at the first blank if absent (defeating the F10.4 implied-decimals rule); parse errors yield 0.
- **Formatted output**: Fortran `Iw`/`Fw.d` editing is replicated, including blank padding, the trailing `.` of `Fw.0`, dropping a leading zero when a field would overflow (`0.74` → `.74`), asterisks on overflow, and `a10`/`a20` blank-padded strings. Console `read` with `f8.4` keeps its implied-decimals and error-branch semantics; output files must not already exist (`status='new'`).
- **Sorting**: the Fortran heapsorts order by key ascending (with payload tie-breaking in `sort2`); equivalent total-order sorts are used, which give identical permutations.
- **Preserved quirks/bugs** (harmless; kept for fidelity and flagged with `NB:` comments):
  - the progress message in the rescaling loop uses the leftover loop variable `ii` rather than `i` (console output only);
  - `neighsim`'s progress message prints the *previous* distance record's names;
  - in `fresca`, `alpha` is updated once more after the converged `phi` is computed, so `Freq_1` uses the post-update value while `Phi_out` reports the pre-update one;
  - the exclusion-file open error path re-reads silently with no message;
  - if `fresca`'s rescaling loop fails to converge within `irepmx` (100) iterations, the reported `Iter` is 101, not 100: Fortran's `DO ir=1,irepmx` loop variable is incremented and tested *before* the loop is abandoned, so it ends one past the limit on normal completion (confirmed against gfortran with an isolated repro of the pattern).
- **Line endings**: outputs use `\n` (Unix convention). The Windows reference files use `\r\n`; strip CR before diffing. Input files may have either ending (CR is stripped on read, emulating Windows text mode).
- **Limits**: the original array bounds are kept (4000 samples, 2000 species, 100 time periods, 2 000 000 observations, 500 000 weights for `frescalo`; 400 000 locations for `sampdist`; 4000 sites, 10 000 species, 5 000 000 records for `neighsim`), including the original error messages and the "press <RETURN> to exit" behaviour.

One deliberate deviation: where the original would read out of bounds or corrupt memory (e.g. `neigh > m` in `sampdist`), the port clamps or fails cleanly instead.