/// [Functions located in the `Faster versions of common Bessel functions.`](<https://docs.scipy.org/doc/scipy/reference/special.html#faster-versions-of-common-bessel-functions>)
mod faster_bessel;

/// Various combinatoric functions for integer-types.
mod combinatorics;

pub use combinatorics::*;
pub use faster_bessel::*;
