use num_traits::real::Real;

/// Modified Bessel function of order 0.
///
/// ## Notes
/// * The range is partitioned into the two intervals [0, 8] and (8, infinity).
/// * [Scipy has this as a
/// ufunc](<https://docs.scipy.org/doc/scipy/reference/special.html#special-functions-scipy-special>),
/// as a supposed wrapper over the Cephes routine.
pub fn i0<F>(x: F) -> F
where
    F: Real,
    f64: From<F>,
{
    return F::from(unsafe { special_fun::unsafe_cephes_double::i0(f64::from(x)) }).unwrap();
}

/// Exponentially scaled modified Bessel function of order 0.
///
/// ## Notes
/// * The range is partitioned into the two intervals [0, 8] and (8, infinity).
/// * [Scipy has this as a
/// ufunc](<https://docs.scipy.org/doc/scipy/reference/special.html#special-functions-scipy-special>),
/// as a supposed wrapper over the Cephes routine.
pub fn i0e<F>(x: F) -> F
where
    F: Real,
    f64: From<F>,
{
    return F::from(unsafe { special_fun::unsafe_cephes_double::i0e(f64::from(x)) }).unwrap();
}
