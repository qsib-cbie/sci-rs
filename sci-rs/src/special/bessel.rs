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
