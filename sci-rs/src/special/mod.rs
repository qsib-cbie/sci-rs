//! Special mathematical functions
//!
//! # Available Functions
//! - Factorial, double factorial, and `k`-factorial
//! - Combinatorics (choice and permutations)

mod combinatorics;
mod factorial;

pub use combinatorics::*;
pub use factorial::*;

/// Adds the [Bessel] trait.
mod bessel;
pub use bessel::Bessel;

// Name is from special/xsf folder, which is Scipy has designated as X special functions (written
// in C++) that are not exposed to Python. We keep these set of functions as being crate internal.
pub(crate) mod xsf;
