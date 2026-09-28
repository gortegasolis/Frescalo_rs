//! Rust port of Mark Hill's FRESCALO suite (Frescalo_1.f, neighsim_1.f,
//! sampdist_1.f, January-June 2011; Hill 2012, *Methods in Ecology and
//! Evolution* 3: 195–205).
//!
//! FRESCALO estimates how the frequency of each species changed over time
//! while allowing for uneven recording effort. The three programs run in
//! sequence:
//!
//! 1. `sampdist` ranks each site's geographic neighbours
//!    ([`neighbourhood::spatial`], [`distance`]).
//! 2. `neighsim` turns spatial ranks plus floristic similarity into
//!    neighbourhood weights ([`neighbourhood::similarity`],
//!    [`neighbourhood::weights`]).
//! 3. `frescalo` pools records over each neighbourhood into local
//!    frequencies, standardises recording effort ([`effort`]) and estimates a
//!    time factor per species and period ([`trend`]). The encounter model that
//!    links effort to the chance of recording a species lives in
//!    [`detection::poisson`].
//!
//! A plain-language guide for ecologists is in the project README
//! ("How it works: a guide for ecologists").
//!
//! Module map:
//!
//! | Module | Contents |
//! |---|---|
//! | [`config`] | Parameters and their defaults (neighbourhood sizes, Φ, R\*, limits) |
//! | [`cli`] | Command-line options and the interactive-prompt fallback |
//! | [`distance`] | Planar and geodesic (WGS84) distances |
//! | [`neighbourhood`] | Spatial ranking, floristic similarity, neighbourhood weights |
//! | [`dataset`] | Reading the `frescalo` inputs (weights, occurrences, exclusions) |
//! | [`detection`] | The Poisson encounter model |
//! | [`effort`] | Local frequencies, effort standardisation (`fresca`), benchmarks, sampling intensity |
//! | [`trend`] | Time factors (`tfcalc`) |
//! | [`output`] | Writers for every output file |
//! | [`fortran`] | Fortran-emulation plumbing (fixed-width text, formatting, sorting, console) |
//!
//! Everything in the default configuration reproduces a native gfortran build
//! of the original programs byte for byte; see the README's fidelity notes.

pub mod cli;
pub mod config;
pub mod dataset;
pub mod detection;
pub mod distance;
pub mod effort;
pub mod fortran;
pub mod neighbourhood;
pub mod output;
pub mod trend;

pub use fortran::*;
