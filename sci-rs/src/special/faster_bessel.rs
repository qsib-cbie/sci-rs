use ndarray::{Array, Dimension};
use num_traits::real::Real;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// All [functions located in the `Faster versions of common Bessel
/// functions.`](<https://docs.scipy.org/doc/scipy/reference/special.html#faster-versions-of-common-bessel-functions>)
pub trait Bessel {
    /// Modified Bessel function of order 0.
    ///
    /// ## Notes
    /// * The range is partitioned into the two intervals [0, 8] and (8, infinity).
    /// * [Scipy has this as a
    ///   ufunc](<https://docs.scipy.org/doc/scipy/reference/special.html#special-functions-scipy-special>),
    ///   as a supposed wrapper over the Cephes routine. We try to define it over reasonable types in
    ///   the impl.
    fn i0(&self) -> Self;

    /// Exponentially scaled modified Bessel function of order 0.
    ///
    /// ## Notes
    /// * The range is partitioned into the two intervals [0, 8] and (8, infinity).
    /// * [Scipy has this as a
    ///   ufunc](<https://docs.scipy.org/doc/scipy/reference/special.html#special-functions-scipy-special>),
    ///   as a supposed wrapper over the Cephes routine. We try to define it over reasonable types in
    ///   the impl.
    fn i0e(&self) -> Self;
}

impl Bessel for f32 {
    // Known to yield wrong result.
    fn i0(&self) -> Self {
        unsafe { special_fun::unsafe_cephes_single::i0f(*self) }
    }

    fn i0e(&self) -> Self {
        unsafe { special_fun::unsafe_cephes_single::i0ef(*self) }
    }
}

impl Bessel for f64 {
    fn i0(&self) -> Self {
        unsafe { special_fun::unsafe_cephes_double::i0(*self) }
    }

    fn i0e(&self) -> Self {
        unsafe { special_fun::unsafe_cephes_double::i0e(*self) }
    }
}

#[cfg(feature = "alloc")]
impl<F> Bessel for Vec<F>
where
    F: Real + super::Bessel,
{
    fn i0(&self) -> Self {
        self.iter().map(|f| f.i0()).collect()
    }

    fn i0e(&self) -> Self {
        self.iter().map(|f| f.i0e()).collect()
    }
}

#[cfg(feature = "alloc")]
impl<F, D> Bessel for Array<F, D>
where
    F: Real + super::Bessel,
    D: Dimension,
{
    fn i0(&self) -> Self {
        self.map(|f| f.i0())
    }

    fn i0e(&self) -> Self {
        self.map(|f| f.i0e())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    #[ignore = "upstream"]
    fn i0single() {
        let inp: f32 = 0.213;
        // upstream (https://www.moshier.net/#Cephes) has to fix
        assert_abs_diff_eq!(1.0113744522192416, inp.i0());
        assert_abs_diff_eq!(0.8173484705849442, inp.i0e());
    }

    #[test]
    fn i0double() {
        let inp: f64 = 0.213;
        assert_abs_diff_eq!(1.0113744522192416, inp.i0());
        assert_abs_diff_eq!(0.8173484705849442, inp.i0e());
    }
}
