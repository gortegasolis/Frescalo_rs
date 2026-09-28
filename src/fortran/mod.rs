//! Fortran-emulation plumbing shared by the three programs: fixed-width text
//! fields, the input-file reader, sorting, Fortran-style formatted output and
//! the interactive console. None of this is part of the ecological method; it
//! exists so that the port reads and writes files exactly like the originals
//! (Frescalo_1.f, neighsim_1.f, sampdist_1.f, January-June 2011).
//!
//! Fidelity conventions:
//! * Fortran `real` is single precision, so all floating-point work uses `f32`.
//! * Fortran `character*10` names are represented as blank-padded `[u8; 10]`
//!   byte arrays; comparisons are byte-wise, exactly like Fortran character
//!   comparison of equal-length strings.
//! * Fortran `character*30` packed records are `[u8; 30]`.
//! * All arrays keep the Fortran 1-based indexing: element 0 is unused.
//! * Formatted output replicates Fortran I/F edit descriptors, including
//!   blank padding, the trailing "." of F w.0, dropping of a leading zero
//!   when the field would otherwise overflow, and asterisks on overflow.

pub mod array;
pub mod console;
pub mod format;
pub mod sort;
pub mod text;

pub use array::*;
pub use console::*;
pub use format::*;
pub use sort::*;
pub use text::*;

/// Buffered output file.
pub type LogWriter = std::io::BufWriter<std::fs::File>;
