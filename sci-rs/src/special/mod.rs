/// Various combinatoric functions for integer-types.
mod combinatorics;
#[cfg(feature = "std")]
mod softmax;

pub use combinatorics::*;
#[cfg(feature = "std")]
pub use softmax::*;
