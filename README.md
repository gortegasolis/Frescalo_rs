# Frescalo_rs

A faithful, behaviour-preserving Rust port of Mark Hill's FRESCALO suite (Hill, 2012, *Methods in Ecology and Evolution* 3: 195–205).

The port is organised as a library with one module per step of the method, plus three small programs that match the original ones. With the default settings it gives exactly the original results. Neighbourhood sizes and the other settings can now be given as command-line options, and site distances can optionally be geodesic (see [How the software is structured](#how-the-software-is-structured)).

Three programs are ported, as three binaries:

| Binary | Original | Purpose |
|----------------------|--------------------------------|------------------|
| `sampdist` | `sampdist_1.f` | Distances between sample locations (planar, or geodesic on WGS84); writes nearest neighbours per location. |
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

### Command-line options

Every setting can also be given as an option (`--name value` or `--name=value`; `--help` lists them):

| Program | Option | Default | Meaning |
|---|---|---|---|
| `sampdist` | `--locations FILE` | *(prompt)* | Locations file, `site x y` |
| | `--output FILE` | *(prompt)* | Neighbourhood-distance output file |
| | `--neighbours N` | 200 | Nearest sites kept per site, the site itself included |
| | `--distance planar\|geodesic` | `planar` | `planar`: x y are easting northing. `geodesic`: x y are longitude latitude in degrees, and distances are WGS84 geodesics in km |
| `neighsim` | `--training FILE` | *(prompt)* | Training-set records, `site species` |
| | `--distances FILE` | *(prompt)* | `sampdist` output |
| | `--similarity-out FILE` | *(prompt)* | Similarity output file |
| | `--weights-out FILE` | *(prompt)* | Weights output file for `frescalo` |
| | `--neighbours K` | 100 | Most similar sites kept per site |
| `frescalo` | `--log FILE` | *(prompt)* | Log file |
| | `--occurrences FILE` | *(prompt)* | Occurrence records, `site species period` |
| | `--weights FILE` | *(prompt)* | `neighsim` weights file |
| | `--exclusions FILE` | none | Species to exclude from the benchmarks, one per line |
| | `--samples-out FILE`, `--frequencies-out FILE`, `--trends-out FILE` | *(prompt)* | The three output files |
| | `--phi PHI` | 0.74 | Target local frequency Φ (0.50–0.95) |
| | `--benchmark-limit R` | 0.2703 | Benchmark limit R\* (0.08–0.5) |
| all | `--no-hold` | | Do not wait for <RETURN> before exiting |

How options and prompts combine:

- **No options:** the program asks every question the original asks, in the same order. A blank answer at the neighbours prompt takes the default (200 for `sampdist`, 100 for `neighsim`), as a blank answer to Φ and R\* already did.
- **Any value option given:** settings with a default (neighbours, distance method, exclusions, Φ, R\*) take it silently. File names that are not given are still asked for.
- The closing "Press <RETURN> to exit" pause happens only if something was typed at a prompt. `--no-hold` turns it off in every case.
- Output files must not already exist, as in the originals. With options, an existing output or a missing input file is a fatal error (exit status 1) instead of a re-prompt.

The `frescalo` log is the same either way: in option mode the prompt texts, file names and parameter values are still written to it. The defaults give exactly the original results, so, for example,

``` sh
sampdist --locations sites.txt --output dist.txt --no-hold
neighsim --training training.txt --distances dist.txt \
         --similarity-out sim.txt --weights-out weights.txt --no-hold
frescalo --log log.txt --occurrences records.txt --weights weights.txt \
         --samples-out samples.txt --frequencies-out freq.txt --trends-out trends.txt --no-hold
```

reproduces Hill's settings (N = 200, K = 100, Φ = 0.74, R\* = 0.2703).

## How it works: a guide for ecologists

FRESCALO answers one question: *has this species become more or less frequent over time, once we allow for the fact that some places and periods were recorded much more thoroughly than others?* Atlas and citizen-science data are full of this kind of bias. A well-visited grid square can list twice as many species as the one next to it for no ecological reason. FRESCALO deals with this by comparing each site with its **neighbourhood**, a set of nearby, ecologically similar sites. It then uses the commonest species in that neighbourhood as a yardstick for how hard each site was searched.

The three programs run in sequence, each feeding the next:

```
site coordinates ──► sampdist ──► nearest neighbours ──┐
                                                       ├──► neighsim ──► neighbourhood weights ──┐
reference species list per site ───────────────────────┘                                        │
                                                                                                 ├──► frescalo ──► effort, rescaled frequencies, trends
occurrence records (site, species, period) ──────────────────────────────────────────────────────┘
```

### Step 1: `sampdist` finds each site's geographic neighbours

**Input:** one line per site, `site easting northing` (projected coordinates, e.g. metres on a national grid).

For every site the program measures the straight-line distance to every other site and keeps the *N* closest, including the site itself at rank 1. Hill (2012) used N = 200 for British hectads, which is the default (`--neighbours` changes it). This is only a shortlist; floristic similarity narrows it down in step 2.

**Sites in longitude/latitude.** Straight-line distance only makes sense on a projected grid. On raw degrees it is wrong in two ways: a degree of longitude shrinks towards the poles (it is 111 km at the equator but only 56 km at 60°), and sites either side of the 180° meridian look half a world apart. With `--distance geodesic`, give `site longitude latitude` in decimal degrees instead. Distances are then the shortest paths over the WGS84 ellipsoid (Karney's algorithm, via `geographiclib-rs`), in kilometres. `sampdist` prints a warning if planar mode is used on coordinates that all look like degrees.

Only the *order* of neighbours is used later on, never the distances themselves. So the choice matters whenever it changes which sites make the shortlist or their order. On the test grid in `tests/fixtures/Samp_lonlat.txt` (100 sites from 35° to 66.5° N, straddling the 180° meridian), planar-on-degrees and geodesic give different 10-site shortlists for 92 of the 100 sites, and on average share only 76% of their members (63% for 25-site shortlists). That grid is an extreme case. But even without the 180° meridian, a planar shortlist on degrees is stretched north–south on the ground, by a factor of about 1.7 at 55° N and 2 at 60° N.

**Output:** `site  neighbour  spatial_rank  distance` (distance in input units for planar mode, km for geodesic mode).

### Step 2: `neighsim` makes the neighbourhoods and their weights

**Inputs:** a *training set* of `site species` records, and the `sampdist` output. The training set is a well-recorded group used only to judge which sites are ecologically alike. It can be the group you are analysing, with all periods pooled, or a better-recorded group.

1. **Floristic similarity.** For each pair of sites it computes the Sørensen index: twice the number of shared species divided by the sum of the two sites' species counts (0 = nothing in common, 1 = identical lists).
2. **Choosing neighbours.** Among each site's geographic shortlist, candidates are ranked from most to least similar, and the top *K* are kept (K = 100 in Hill 2012, the default; `--neighbours` changes it). The result is a neighbourhood of sites that are both close by and have a similar flora or fauna.
3. **Weighting.** Close, similar neighbours should count for more than marginal ones. Each neighbour gets a weight
   `w = (1 − r_sim²)⁴ × (1 − r_dist²)⁴`,
   where `r_sim` and `r_dist` are its similarity rank and distance rank scaled from 0 (the site itself) to about 1 (the edge of the neighbourhood). The site itself gets weight 1. Weights fall off smoothly with rank, with no hard cut-off, and weights below 0.00005 are dropped.

**Outputs:** a similarity file (`site  neighbour  similarity_rank  Sørensen`) and the **weights file** used by `frescalo` (`site  neighbour  weight  similarity_weight  distance_weight  K  N`).

### Step 3: `frescalo` corrects for recording effort and estimates trends

**Inputs:** the occurrence records (`site species period`), the weights file, and optionally a list of species that should not be used as benchmarks. There are also two parameters, **Φ** (phi, default 0.74, `--phi`) and the **benchmark limit R\*** (default 0.27, `--benchmark-limit`). The program works in four stages.

**(a) Local frequencies.** Records from all periods are pooled. For each site and species, the program computes the weighted proportion of the site's neighbourhood where that species was found. This is the species' *local frequency*, *f*. A species found in every neighbouring square has *f* ≈ 1, and one found in a single, low-weight neighbour has *f* close to 0. It is a smoothed estimate of how likely the species is to occur at that site.

**(b) Standardising effort per neighbourhood: `fresca`.** Neighbourhoods differ in both richness and recording effort. `fresca` separates the two, one site at a time. It summarises the neighbourhood by Φ, the *frequency-weighted mean frequency* (Σf² / Σf). Heavily recorded neighbourhoods have high Φ, because their common species really do turn up nearly everywhere. Under-recorded neighbourhoods have low Φ.

`fresca` then looks for the **sampling-effort multiplier α** that would bring the neighbourhood's Φ to the common target value. It does this by "turning up the effort" on every species at once:

   `f_rescaled = 1 − (1 − f)^α`

   This formula comes from a Poisson model of recording. If a recorder meets the species at random, at an average rate λ during a visit, the chance of recording it at least once is `1 − e^(−λ)`. Multiplying the effort by α multiplies λ by α. The model is in its own module (`detection::poisson`), shared by `fresca` and `tfcalc`.

It tries a value of α, checks the resulting Φ, adjusts, and repeats until Φ is within 0.0003 of the target (at most 100 tries). The rescaled frequencies are the chance of finding each species at that site *if it had been recorded to the common standard*. With every site on the same footing, sites can be compared fairly.

   - α > 1: the neighbourhood is under-recorded compared with the standard, so its frequencies are scaled up.
   - α < 1: it is recorded more thoroughly than the standard, so they are scaled down.

`fresca` also ranks the species at each site from most to least frequent. The most frequent species, those in roughly the top 27% (R\*) of the expected species count, become that site's **benchmark species**. These are species that *should* be found by anyone recording there. Species on the exclusion list are kept as benchmarks but count almost nothing (weight 0.001). Use the list for common species that are easily missed or have changed genuinely.

**(c) Recording effort per site and period.** For each site and time period, the program works out what fraction of that site's benchmark species were actually recorded in that period. This is the **sampling intensity** of the site in that period. A value of 0.9 means a thorough survey. A value of 0.1 means the site was barely visited, and a species' absence tells us little.

**(d) Time factors: `tfcalc`.** For each species and period, `tfcalc` combines each site's local frequency with that period's sampling intensity. This gives the probability that the species *would* have been recorded there if nothing had changed. Adding these probabilities over all sites gives the *expected* number of records. `tfcalc` then finds the **time factor** that makes the expected number match the observed number, again by repeated adjustment. Sites with very low sampling intensity (below about 0.1) are down-weighted, because their records carry little information.

   - Time factor 1: the species occurred in that period about as often as its pooled local frequency predicts.
   - Above 1: it was more frequent than that, relative to effort.
   - Below 1: it was less frequent.

   Time factors are relative. What matters ecologically is how a species' time factor *changes across periods*, not its absolute value. A **standard error** comes from asking how far the time factor would move if the observed count were one standard deviation higher.

When the run finishes, the log reports the 98.5th percentile of the sites' raw Φ values next to your target Φ. If the target is lower than this, the program prints a warning. A target that low means most of the best-recorded neighbourhoods would have to be scaled *down* (α < 1). Hill (2012) recommends a target near the top of the observed range, and the 0.74 default was chosen for British bryophytes.

### Reading the `frescalo` output files

**Sample statistics** (one line per site; the output of `fresca`):

| Column | Meaning |
|---|---|
| `Location`, `Loc_no` | Site code and its serial number |
| `No_spp` | Number of species actually recorded at the site (all periods) |
| `Phi_in` | Neighbourhood Φ before rescaling; low values mean under-recorded |
| `Alpha` | Sampling-effort multiplier α (capped at 999.99 in the output); > 1 means under-recorded |
| `Wgt_n2` | Effective number of neighbours (few neighbours means less reliable frequencies) |
| `Phi_out` | Φ after rescaling; should equal the target |
| `Spnum_in` | Expected number of species at the site at its actual effort (sum of local frequencies) |
| `Spnum_out` | Expected number of species at standard effort, a bias-corrected richness estimate |
| `Iter` | Iterations needed to find α; 101 means it failed to converge |

**Rescaled frequencies** (one line per site × species):

| Column | Meaning |
|---|---|
| `Location`, `Species` | Site and species |
| `Pres` | 1 if the species was recorded at the site itself, 0 otherwise |
| `Freq__` | Raw local frequency *f* |
| `Freq_1` | Rescaled frequency: chance of recording the species here at standard effort |
| `SD_Frq1` | Standard error of `Freq_1` (larger when the neighbourhood is small) |
| `Rank`, `Rank_1` | Rank of the species at this site, and rank ÷ `Spnum_out`; species with `Rank_1` below R\* are benchmarks |

**Trends** (one line per species × period):

| Column | Meaning |
|---|---|
| `Species`, `Time` | Species and time period |
| `TFactor`, `St_Dev` | Time factor and its standard error |
| `_Count` | Number of sites where the species was recorded in that period |
| `___spt` | The same count, with poorly sampled sites (intensity below about 0.1) down-weighted |
| `___est` | The expected weighted count at the fitted time factor (≈ `___spt` when converged) |
| `N>0.00` | Sites where the species had any chance of being recorded in that period |
| `N>0.98` | Sites where the chance was essentially certain (capped at 0.98) |

Periods with no records for a species are written as rows of zeros, so every species has one row per period.

## How the software is structured

The original FRESCALO is three self-contained Fortran programs. Each is one long main routine that reads its answers from the keyboard, with the numbers of neighbours and other settings fixed in the code. The first version of this port copied that shape closely. The current version keeps the same calculations but arranges them differently:

- **One library, three thin programs.** All of the method lives in a library (`src/lib.rs` and its modules). The programs in `src/bin/` (`sampdist.rs`, `neighsim.rs`, `frescalo.rs`) only read options or prompts, open files, call the library in the right order, and report errors.
- **One module per step of the method.** Each step described in the guide above has its own module, so the code for a step can be read, tested and changed without wading through the rest.
- **Settings in one place.** Every parameter and its default lives in `config`. The neighbourhood sizes that were fixed at 200 and 100 are now options that default to those values.
- **The science kept apart from the Fortran emulation.** Everything needed to reproduce the Fortran output exactly (fixed-width names, number formatting, sort order, console behaviour) sits in `fortran`. The method modules call it but do not reimplement it.

```
                 src/bin/sampdist.rs   src/bin/neighsim.rs   src/bin/frescalo.rs
programs         options, prompts and files (cli, config)
                        │                     │                     │
                        ▼                     ▼                     ▼
method           distance            neighbourhood::        dataset → effort → trend
                 neighbourhood::       similarity                     │        │
                   spatial             weights                        └─ detection::poisson
                        │                     │                     │
                        ▼                     ▼                     ▼
output           output (one writer per file layout)
plumbing         fortran (text fields, formatting, sorting, console)
```

### Where each step lives

| Step of the method (see the guide above) | Module | Main functions |
|---|---|---|
| Parameters and their defaults | `config` | `SampdistParams`, `NeighsimParams`, `FrescaloParams`, `limits` |
| Options and prompts | `cli` | `Args`, `Session` |
| 1. Distances between sites | `distance` | `Metric` (planar or geodesic), `DistanceEngine` |
| 1. Each site's nearest sites | `neighbourhood::spatial` | `read_locations`, `rank_by_distance` |
| 2. Floristic similarity | `neighbourhood::similarity` | `read_inputs`, `shared_species`, `sorensen_row` |
| 2. Neighbourhood weights | `neighbourhood::weights` | `mark_spatial_shortlist`, `write_neighbourhoods`, `rank_kernel` |
| 3. Reading records and weights | `dataset` | `Dataset`, `read_exclusions` |
| 3a. Local frequencies | `effort` | `weight_totals`, `local_frequencies` |
| 3b. Standardising effort | `effort` | `fresca`, `standardise_sites` |
| 3c. Sampling intensity | `effort` | `sampling_intensity` |
| 3d. Time factors | `trend`, `trend::hill` | `write_trends`, `tfcalc` |
| Chance of recording a species, given effort | `detection::poisson` | `intensity`, `prob_from_intensity`, `rescale`, `alpha_step` |
| Output files | `output` | one `write_*` function per file layout |
| Fortran emulation | `fortran` | text fields, `ffmt`/`ifmt`, sorting, console input |

### Why the Poisson model has its own module

The assumption that recording follows a Poisson process (the chance of recording a species is `1 − e^(−λ)`, and extra effort multiplies λ) is the core statistical assumption of FRESCALO. In the Fortran it is written out inline in five places across `fresca` and `tfcalc`. It is now defined once, in `detection::poisson`, with a plain-language explanation. `effort` and `trend` call it rather than repeating the formulas. This makes the assumption easy to find and check, and means an alternative detection model could be added next to it without rewriting the effort and trend code.

### What stays the same

- **Results.** Run with the default settings, the programs produce output byte-identical to a gfortran build of the original Fortran. The restructuring moved code but did not change the order of any calculation. `verify_against_fortran.sh` and the tests in `tests/golden.rs` check this.
- **The interactive programs.** Run without options, each program asks the original questions in the original order. Existing scripts that pipe answers into the programs keep working.
- **The file formats.** Inputs and outputs are unchanged. The only addition is geodesic mode in `sampdist`, which reads longitude and latitude and writes distances in km.

### Extending the software

- **Change a default:** edit `config.rs`. The defaults are named constants, and the help text of the matching option in `src/bin/` should say the same.
- **Add a distance method:** add a variant to `Metric` in `distance.rs`, its calculation to `DistanceEngine::distance`, and its name to `Metric::from_option` and the `--distance` help text in `src/bin/sampdist.rs`. Nothing downstream changes, because `neighsim` uses only the order of neighbours.
- **Add a detection model:** add a module next to `detection/poisson.rs`, and change the calls to `poisson::` in `effort.rs` and `trend/hill.rs`. The Poisson functions must keep being used unchanged when the default is selected, or the outputs will no longer match the Fortran.
- **Check a change:** `cargo test --release` runs the unit tests in each module and the golden-output tests. Any change that alters results in the default configuration makes the golden tests fail.

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
./frescalo.sh verify    # run verify_against_fortran.sh (needs a folder named Original_frescalo in the root, containing the original fortran code files)
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

`cargo test` needs no Fortran compiler. `tests/golden.rs` runs the pipeline on a small synthetic data set (`tests/fixtures/`, made by `generate.py`) and compares every output with files produced by the gfortran build. It covers prompted runs, option-driven runs, runs mixing the two, the defaults for a blank answer, and geodesic mode.

### Relationship to the Windows reference outputs

The Rust port reproduces results from the original fortran code exactly except for last-digit rounding in a small fraction of lines (`dist.txt` 0/80800, `samples.txt` 0/405, `trends.txt` 3/6521, `frequencies.txt` 583/185803, `sim.txt` 16/40400, `weights.txt` 12/35549). These differences are platform artifacts of the old Windows compiler: it evaluated intermediate expressions in x87 80-bit extended precision and rounded exact ties half-away-from-zero, whereas modern x86-64 (both gfortran and this port) uses strict IEEE single precision and round-half-to-even. A native gfortran build of the original sources produces output identical to this port, not to the Windows reference files.

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
- **Limits**: the original array bounds are kept (4000 samples, 2000 species, 100 time periods, 2 000 000 observations, 500 000 weights for `frescalo`; 400 000 locations for `sampdist`; 4000 sites, 10 000 species, 5 000 000 records for `neighsim`), including the original error messages. After such an error the program pauses for <RETURN> if it was run interactively, then exits with status 1.

Deliberate deviations:

- Where the original would read out of bounds or corrupt memory (e.g. `neigh > m` in `sampdist`), the port clamps or fails cleanly instead.
- The command-line options, and the defaults for a blank neighbours answer, are additions. The original's list-directed read would skip a blank line and wait for another. None of these change the results obtained with the original inputs.
- **Geodesic mode** (`--distance geodesic`) has no Fortran counterpart. It computes in double precision (`f64`), so it is not bit-for-bit comparable with anything. The distance column is still written with the original `f6.0` format (whole kilometres; the longest possible geodesic, about 20 004 km, fits). Planar mode, the default, keeps the original single-precision arithmetic exactly.