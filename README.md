# Frescalo_rs

A faithful, behaviour-preserving Rust port of Mark Hill's FRESCALO suite (Hill, 2011, *Methods in Ecology and Evolution* 2: 502–512).

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

For every site the program measures the straight-line distance to every other site and keeps the *N* closest, including the site itself at rank 1. Hill (2011) used N = 200 for British hectads. This is only a shortlist; floristic similarity narrows it down in step 2.

**Output:** `site  neighbour  spatial_rank  distance`.

### Step 2: `neighsim` makes the neighbourhoods and their weights

**Inputs:** a *training set* of `site species` records, and the `sampdist` output. The training set is a well-recorded group used only to judge which sites are ecologically alike. It can be the group you are analysing, with all periods pooled, or a better-recorded group.

1. **Floristic similarity.** For each pair of sites it computes the Sørensen index: twice the number of shared species divided by the sum of the two sites' species counts (0 = nothing in common, 1 = identical lists).
2. **Choosing neighbours.** Among each site's geographic shortlist, candidates are ranked from most to least similar, and the top *K* are kept (K = 100 in Hill 2011). The result is a neighbourhood of sites that are both close by and have a similar flora or fauna.
3. **Weighting.** Close, similar neighbours should count for more than marginal ones. Each neighbour gets a weight
   `w = (1 − r_sim²)⁴ × (1 − r_dist²)⁴`,
   where `r_sim` and `r_dist` are its similarity rank and distance rank scaled from 0 (the site itself) to about 1 (the edge of the neighbourhood). The site itself gets weight 1. Weights fall off smoothly with rank, with no hard cut-off, and weights below 0.00005 are dropped.

**Outputs:** a similarity file (`site  neighbour  similarity_rank  Sørensen`) and the **weights file** used by `frescalo` (`site  neighbour  weight  similarity_weight  distance_weight  K  N`).

### Step 3: `frescalo` corrects for recording effort and estimates trends

**Inputs:** the occurrence records (`site species period`), the weights file, and optionally a list of species that should not be used as benchmarks. There are also two parameters, **Φ** (phi, default 0.74) and the **benchmark limit R\*** (default 0.27). The program works in four stages.

**(a) Local frequencies.** Records from all periods are pooled. For each site and species, the program computes the weighted proportion of the site's neighbourhood where that species was found. This is the species' *local frequency*, *f*. A species found in every neighbouring square has *f* ≈ 1, and one found in a single, low-weight neighbour has *f* close to 0. It is a smoothed estimate of how likely the species is to occur at that site.

**(b) Standardising effort per neighbourhood: `fresca`.** Neighbourhoods differ in both richness and recording effort. `fresca` separates the two, one site at a time. It summarises the neighbourhood by Φ, the *frequency-weighted mean frequency* (Σf² / Σf). Heavily recorded neighbourhoods have high Φ, because their common species really do turn up nearly everywhere. Under-recorded neighbourhoods have low Φ.

`fresca` then looks for the **sampling-effort multiplier α** that would bring the neighbourhood's Φ to the common target value. It does this by "turning up the effort" on every species at once:

   `f_rescaled = 1 − (1 − f)^α`

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

When the run finishes, the log reports the 98.5th percentile of the sites' raw Φ values next to your target Φ. If the target is lower than this, the program prints a warning. A target that low means most of the best-recorded neighbourhoods would have to be scaled *down* (α < 1). Hill (2011) recommends a target near the top of the observed range, and the 0.74 default was chosen for British bryophytes.

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
- **Limits**: the original array bounds are kept (4000 samples, 2000 species, 100 time periods, 2 000 000 observations, 500 000 weights for `frescalo`; 400 000 locations for `sampdist`; 4000 sites, 10 000 species, 5 000 000 records for `neighsim`), including the original error messages and the "press <RETURN> to exit" behaviour.

One deliberate deviation: where the original would read out of bounds or corrupt memory (e.g. `neigh > m` in `sampdist`), the port clamps or fails cleanly instead.