use num_traits::real::Real;

/// This function evaluates the series y = Sum{i=0..N} coef[i] T_i(x/2) of Chebyshev polynomials Ti
/// at argument x/2. Returns 0 for empty coef.
///
/// Coefficients are stored in the reverse order, i.e.: the zero order term is the last in the
/// array.
pub(crate) fn chbevl<T>(x: T, coef: &[T]) -> T
where
    T: Real + core::ops::Mul,
{
    // Summation over the empty set is defined to be zero.
    if coef.is_empty() {
        return T::zero();
    }

    let (b0, _, b2): (T, T, T) = coef.iter().fold(
        (
            // Safety:: 0 len is checked above.
            *unsafe { coef.first().unwrap_unchecked() },
            T::zero(),
            T::zero(),
        ),
        |acc, &e| (x * acc.0 - acc.1 + e, acc.0, acc.1),
    );

    (b0 - b2) / unsafe { T::from(2.).unwrap_unchecked() }
}
