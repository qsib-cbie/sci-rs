use nalgebra::min;
use num_traits::{FromPrimitive, PrimInt};

/// The number of combinations of `n` taken `k` at a time.
///
/// This is also known as $n$ choose $k$ and is generally given by the formula
/// $$
/// \begin{pmatrix}
/// n \\\\ k
/// \end{pmatrix} = \frac{n!}{k!(n-k)!}
/// $$
///
/// # Examples
/// ```
/// use sci_rs::special::comb;
/// assert_eq!(comb(5, 2), 10);
/// ```
///
/// # Notes
/// When `n < 0` or `k < 0` or `n < k`, then `0` is returned.
pub fn comb<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if k > n || n < Int::zero() || k < Int::zero() {
        return Int::zero();
    }
    let m = n + Int::one();
    let n_terms = (min(k, n - k) + Int::one()).to_usize().unwrap();
    (1..n_terms).fold(Int::one(), |result, i| {
        result * (m - Int::from_usize(i).unwrap()) / Int::from_usize(i).unwrap()
    })
}

/// Number of combinations with repetition.
///
/// This is also known as a `k`-combination with repetition or `k`-multicombinations. For a more
/// detailed explanation, see the [wiki] page.
///
/// # Examples
/// ```
/// use sci_rs::special::comb_rep;
/// assert_eq!(comb_rep(5, 2), 15);
/// assert_eq!(comb_rep(10, 3), 220);
/// ```
///
/// # References
/// - [wiki]
///
/// [wiki]: https://en.wikipedia.org/wiki/Combination#Number_of_combinations_with_repetition
pub fn comb_rep<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    comb(n + k - Int::one(), k)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::fmt;

    fn check_values<T>(x: T, ref_values: &[T], func: fn(T, T) -> T)
    where
        T: PrimInt + FromPrimitive + fmt::Debug,
    {
        for (i, &val) in ref_values.iter().enumerate() {
            let i = T::from_usize(i).unwrap();
            assert_eq!(func(x, i), val);
        }
    }

    #[test]
    fn choose() {
        assert_eq!(comb(3_u8, 1), 3);
        assert_eq!(comb(3_u8, 2), 3);
        assert_eq!(comb(3_u8, 3), 1);

        const REF_VALUES_5: [i32; 7] = [1, 5, 10, 10, 5, 1, 0];
        const REF_VALUES_10: [i32; 12] = [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1, 0];
        const REF_VALUES_15: [i32; 16] = [
            1, 15, 105, 455, 1365, 3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15, 1,
        ];
        check_values(5, &REF_VALUES_5, comb);
        check_values(10, &REF_VALUES_10, comb);
        check_values(15, &REF_VALUES_15, comb);
    }

    #[test]
    fn choose_negatives() {
        for i in 0..4 {
            assert_eq!(comb(-4, i), 0);
            assert_eq!(comb(-3, i), 0);
            assert_eq!(comb(-3241, i), 0);
        }
        for i in -4..0 {
            assert_eq!(comb(4, i), 0);
            assert_eq!(comb(2, i), 0);
            assert_eq!(comb(2341, i), 0);
            assert_eq!(comb(-2, i), 0);
            assert_eq!(comb(-4, i), 0);
            assert_eq!(comb(-5, i), 0);
            assert_eq!(comb(-3241, i), 0);
        }
    }

    #[test]
    fn zero_choose_zero() {
        assert_eq!(comb(0, 0), 1);
    }
}
