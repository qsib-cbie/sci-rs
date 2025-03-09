/// Various combinatoric functions for integer-types.
mod combinatorics;

pub use combinatorics::*;

// Name is from special/xsf folder, which is Scipy has designated as X special functions (written
// in C++) that are not exposed to Python. We keep these set of functions as being crate internal.
pub(crate) mod xsf;
