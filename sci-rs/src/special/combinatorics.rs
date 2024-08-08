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
/// # Notes
/// When `n < 0` or `k < 0` or `n < k`, then `0` is returned.
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
        for n in -10..-1 {
            for m in -5..5 {
                assert_eq!(comb(n, m), 0);
                assert_eq!(comb(m, n), 0);
            }
        }
    }

    #[test]
    fn zero_choose_zero() {
        assert_eq!(comb(0, 0), 1);
    }

    #[test]
    fn choose_replacement() {
        const REF_VALUES_5: [i32; 10] = [1, 5, 15, 35, 70, 126, 210, 330, 495, 715];
        const REF_VALUES_7: [i32; 10] = [1, 7, 28, 84, 210, 462, 924, 1716, 3003, 5005];
        const REF_VALUES_10: [i32; 15] = [
            1, 10, 55, 220, 715, 2002, 5005, 11440, 24310, 48620, 92378, 167960, 293930, 497420,
            817190,
        ];
        check_values(5, &REF_VALUES_5, comb_rep);
        check_values(7, &REF_VALUES_7, comb_rep);
        check_values(10, &REF_VALUES_10, comb_rep);
    }

    #[test]
    fn choose_replacement_negatives() {
        for n in -10..-1 {
            for m in -5..5 {
                assert_eq!(comb_rep(n, m), 0);
                assert_eq!(comb_rep(m, n), 0);
            }
        }
    }

    #[test]
    fn choose_zero_replacement() {
        for i in 0..1 {
            for j in 0..1 {
                assert_eq!(comb_rep(i, j), i);
            }
        }
    }
}
