use ndarray::ArrayD;
use num_traits::real::Real;

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
    ///
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

#[cfg(feature = "std")]
impl<T> Bessel for Vec<T>
where
    T: Bessel,
{
    fn i0(&self) -> Self {
        self.iter().map(Bessel::i0).collect()
    }

    fn i0e(&self) -> Self {
        self.iter().map(Bessel::i0e).collect()
    }
}

#[cfg(feature = "std")]
impl<T> Bessel for ArrayD<T>
where
    T: Bessel,
{
    fn i0(&self) -> Self {
        self.map(Bessel::i0)
    }

    fn i0e(&self) -> Self {
        self.map(Bessel::i0e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[cfg(feature = "std")]
    #[test]
    fn i0_vec_f64() {
        let inp = vec![0., 1., 0.213, 5., 30.546];
        let result = vec![
            1.,
            1.2660658777520082,
            1.0113744522192416,
            27.239871823604442,
            1337209608661.4026,
        ];
        assert_eq!(result.len(), inp.i0().len());
        for (&r, e) in result.iter().zip(inp.i0()) {
            assert_relative_eq!(r, e, epsilon = 1e-6);
        }
    }
}
