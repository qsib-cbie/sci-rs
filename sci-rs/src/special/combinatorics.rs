use nalgebra::min;
use num_traits::{FromPrimitive, PrimInt};

/// Various combinatorics functions for integer types.
pub trait Combinatoric {
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
    /// use sci_rs::special::Combinatoric;
    /// assert_eq!(5_i32.comb(2), 10);
    /// ```
    ///
    /// # Notes
    /// When `n < 0` or `k < 0` or `n < k`, then `0` is returned.
    fn comb(self, k: Self) -> Self;

    /// Number of combinations with repetition.
    ///
    /// The number of combinations of `n` taken `k` at a time with repetition. This is also known as a
    /// `k`-combination with repetition or `k`-multicombinations. For a more detailed explanation, see
    /// the [wiki] page.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Combinatoric;
    /// assert_eq!(5_i32.comb_rep(2), 15);
    /// assert_eq!(10_i32.comb_rep(3), 220);
    /// ```
    ///
    /// # Notes
    /// When `n < 0` or `k < 0` or `n < k`, then `0` is returned.
    ///
    /// # References
    /// - [Wikipedia][wiki]
    ///
    /// [wiki]: https://en.wikipedia.org/wiki/Combination#Number_of_combinations_with_repetition
    fn comb_rep(self, k: Self) -> Self;

    /// Number of permutations of `n` things taken `k` at a time.
    ///
    /// Also known as the `k`-permutation of `n`.
    /// $$
    /// \text{Perm}(n, k) = \frac{n!}{(n-k)!}
    /// $$
    /// # Examples
    /// ```
    /// use sci_rs::special::Combinatoric;
    /// assert_eq!(5.perm(5), 120); // should be 5!
    /// assert_eq!(5.perm(0), 1);
    /// assert_eq!(6.perm(3), 6*5*4);
    /// ```
    ///
    /// # Notes
    /// When `n<0` or `k<0`, the `0` is returned.
    fn perm(self, k: Self) -> Self;

    /// Stirling number of the second kind.
    ///
    /// These count the number of ways to partition a set of `n` elements into `k` non-empty
    /// subsets. These are often called `n` subset `k` and denoted as either
    /// $$
    /// S(n,k)
    /// $$
    /// or
    /// $$
    /// \begin{Bmatrix}
    ///     n \\ k
    /// \end{Bmatrix}
    /// $$
    /// See the [wiki] page for more details.
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Combinatoric;
    /// assert_eq!(3.stirling2(2), 3);
    /// assert_eq!(0.stirling2(0), 1);
    /// assert_eq!(4.stirling2(3), 6);
    /// ```
    ///
    /// # References
    /// - [Wikipedia][wiki]
    ///
    /// [wiki]: https://en.wikipedia.org/wiki/Stirling_numbers_of_the_second_kind
    fn stirling2(self, k: Self) -> Self;
}

macro_rules! combinatoric_primint_impl {
    ($($T: ty)*) => ($(
        impl Combinatoric for $T {
            #[inline(always)]
            fn comb(self, k: Self) -> Self {
                primint_comb(self, k)
            }

            #[inline(always)]
            fn comb_rep(self, k: Self) -> Self {
                primint_comb_rep(self, k)
            }

            #[inline(always)]
            fn perm(self, k: Self) -> Self {
                primint_perm(self, k)
            }

            #[inline(always)]
            fn stirling2(self, k: Self) -> Self {
                primint_stirling2(self, k)
            }

        }
    )*)
}

combinatoric_primint_impl! {u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize}

fn primint_comb<Int>(n: Int, k: Int) -> Int
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

fn primint_comb_rep<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if n + k == Int::zero() {
        return Int::zero();
    }
    primint_comb(n + k - Int::one(), k)
}

fn primint_perm<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if k > n || n < Int::zero() || k < Int::zero() {
        return Int::zero();
    }

    let start = (n - k + Int::one()).to_usize().unwrap();
    let end = (n + Int::one()).to_usize().unwrap();
    (start..end).fold(Int::one(), |result, val| {
        result * Int::from_usize(val).unwrap()
    })
}

fn primint_stirling2<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if n < Int::zero() || k < Int::zero() {
        return Int::zero();
    }
    if k > n {
        return Int::zero();
    }

    if n == k {
        return Int::one();
    }

    if k == Int::zero() || n == Int::zero() {
        return Int::zero();
    }

    k * primint_stirling2(n - Int::one(), k) + primint_stirling2(n - Int::one(), k - Int::one())
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::fmt;

    fn check_values<K, T>(ref_values: &[[K;10]], func: fn(T, T) -> T)
    where
        K: PrimInt + FromPrimitive, 
        T: PrimInt + FromPrimitive + fmt::Debug,
    {
        for (n, &elements) in ref_values.iter().enumerate() {
            for (k, &val) in elements.iter().enumerate() {
                let n = T::from_usize(n).unwrap();
                let k = T::from_usize(k).unwrap();
                let val = T::from(val).unwrap();
                assert_eq!(func(n, k), val);
            }
        }
    }

    #[test]
    fn comb() {
        // Generated from scipy
        let ref_values = [
           [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 2, 1, 0, 0, 0, 0, 0, 0, 0],
           [1, 3, 3, 1, 0, 0, 0, 0, 0, 0],
           [1, 4, 6, 4, 1, 0, 0, 0, 0, 0],
           [1, 5, 10, 10, 5, 1, 0, 0, 0, 0],
           [1, 6, 15, 20, 15, 6, 1, 0, 0, 0],
           [1, 7, 21, 35, 35, 21, 7, 1, 0, 0],
           [1, 8, 28, 56, 70, 56, 28, 8, 1, 0],
           [1, 9, 36, 84, 126, 126, 84, 36, 9, 1],
        ];
        check_values(&ref_values, u32::comb);
        check_values(&ref_values, i32::comb);
    }

    #[test]
    fn comb_negatives() {
        for n in -10..0 {
            for m in -5..5 {
                assert_eq!(n.comb(m), 0);
                assert_eq!(m.comb(n), 0);
            }
        }
    }

    #[test]
    fn comb_rep() {
        // Generated from scipy
        let ref_values = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            [1, 3, 6, 10, 15, 21, 28, 36, 45, 55],
            [1, 4, 10, 20, 35, 56, 84, 120, 165, 220],
            [1, 5, 15, 35, 70, 126, 210, 330, 495, 715],
            [1, 6, 21, 56, 126, 252, 462, 792, 1287, 2002],
            [1, 7, 28, 84, 210, 462, 924, 1716, 3003, 5005],
            [1, 8, 36, 120, 330, 792, 1716, 3432, 6435, 11440],
            [1, 9, 45, 165, 495, 1287, 3003, 6435, 12870, 24310],
        ];
        check_values(&ref_values, u32::comb_rep);
        check_values(&ref_values, i32::comb_rep);
    }

    #[test]
    fn comb_rep_negatives() {
        for n in -4..0 {
            for m in -5..5 {
                assert_eq!(n.comb_rep(m), 0);
                assert_eq!(m.comb_rep(n), 0);
            }
        }
    }

    #[test]
    fn perm() {
        // Generated from scipy
        let ref_values = [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 2, 2, 0, 0, 0, 0, 0, 0, 0],
            [1, 3, 6, 6, 0, 0, 0, 0, 0, 0],
            [1, 4, 12, 24, 24, 0, 0, 0, 0, 0],
            [1, 5, 20, 60, 120, 120, 0, 0, 0, 0],
            [1, 6, 30, 120, 360, 720, 720, 0, 0, 0],
            [1, 7, 42, 210, 840, 2520, 5040, 5040, 0, 0],
            [1, 8, 56, 336, 1680, 6720, 20160, 40320, 40320, 0],
            [1, 9, 72, 504, 3024, 15120, 60480, 181440, 362880, 362880],
        ];
        check_values(&ref_values, u32::perm);
        check_values(&ref_values, i32::perm);
    }

    #[test]
    fn perm_negative() {
        for i in -4..0 {
            for j in -5..5 {
                assert_eq!(i.perm(j), 0);
                assert_eq!(j.perm(i), 0);
            }
        }
    }

    #[test]
    fn stirling2() {
        // Generated from scipy
        let ref_values = [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 3, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 7, 6, 1, 0, 0, 0, 0, 0],
            [0, 1, 15, 25, 10, 1, 0, 0, 0, 0],
            [0, 1, 31, 90, 65, 15, 1, 0, 0, 0],
            [0, 1, 63, 301, 350, 140, 21, 1, 0, 0],
            [0, 1, 127, 966, 1701, 1050, 266, 28, 1, 0],
            [0, 1, 255, 3025, 7770, 6951, 2646, 462, 36, 1],
        ];
        check_values(&ref_values, u32::stirling2);
        check_values(&ref_values, i32::stirling2);
    }

    #[test]
    fn stirling2_negative() {
        for i in -4..0 {
            for j in -5..5 {
                assert_eq!(i.stirling2(j), 0);
                assert_eq!(j.stirling2(i), 0);
            }
        }
    }
}
