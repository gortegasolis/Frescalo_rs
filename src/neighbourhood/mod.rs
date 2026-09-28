//! Stage 1: building each site's neighbourhood.
//!
//! [`spatial`] (`sampdist`) ranks every site's geographic neighbours;
//! [`similarity`] and [`weights`] (`neighsim`) keep the floristically most
//! similar of those and weight them by both similarity rank and distance
//! rank. The resulting weights file is what `frescalo` pools records over.

pub mod similarity;
pub mod spatial;
pub mod weights;
