//! Special mathematical functions
//!
//! # Available Functions
//! - Factorial, double factorial, and `k`-factorial
//! - Combinatorics (choice and permutations)

mod combinatorics;
mod factorial;
#[cfg(feature = "std")]
mod softmax;

pub use combinatorics::*;
pub use factorial::*;
#[cfg(feature = "std")]
pub use softmax::*;

/// Adds the [Bessel] trait.
mod bessel;
pub use bessel::Bessel;

// Name is from special/xsf folder, which is Scipy has designated as X special functions (written
// in C++) that are not exposed to Python. We keep these set of functions as being crate internal.
pub(crate) mod xsf;
