//! Detection models: how recording effort turns into the chance of recording
//! a species.
//!
//! FRESCALO uses one model, the Poisson encounter model in [`poisson`]. It is
//! kept in its own module so that the effort standardisation
//! ([`crate::effort`]) and the time factors ([`crate::trend`]) share a single,
//! documented definition, and so that an alternative model could be added
//! later without touching those callers.

pub mod poisson;
