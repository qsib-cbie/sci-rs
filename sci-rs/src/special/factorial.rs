use num_traits::{FromPrimitive, PrimInt};

/// Factorial and related functions
///
/// Note that for primitive integer types, it is known the maximum value one can calculate
///
/// | Type | n! | n!! |
/// |------|----|-----|
/// | u8   | 5  | 7   |
/// | u16  | 8  | 12  |
/// | u32  | 12 | 20  |
/// | u64  | 20 | 33  |
/// | u128 | 35 | 57  |
/// | i8   | 5  | 7   |
/// | i16  | 7  | 11  |
/// | i32  | 12 | 19  |
/// | i64  | 20 | 33  |
/// | i128 | 34 | 57  |
pub trait Factorial {
    /// The factorial functions is defined as the product of all positive integers less than or
    /// equal to `n`.
    ///
    /// $$
    /// n! = n(n-1)(n-2)\times\ldots\times1
    /// $$
    ///
    /// Additionally, we have that `0! = 1`.
    ///
    /// # Examples
    /// ## For `usize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0.factorial(), 1);
    /// assert_eq!(3.factorial(), 6);
    /// assert_eq!(5.factorial(), 120);
    /// ```
    ///
    /// ## For `isize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0_isize.factorial(), 1);
    /// assert_eq!(3_isize.factorial(), 6);
    /// assert_eq!(5_isize.factorial(), 120);
    /// ```
    ///
    /// # Notes
    /// If `n < 0`, returns `0`.
    fn factorial(self) -> Self;

    /// The double factorial `n!!` is defined as the product of all
    /// positive integers up to `n` that have the same parity as `n`.
    /// $$
    /// n!! = n (n-2) (n-4) \times \ldots
    /// $$
    /// In the case that `n` is even (`n = 2k`), then
    /// $$
    /// (2k)!! = (2k)(2k-2)(2k-4)\ldots(4)(2)
    /// $$
    /// In the case that `n` is odd (`n=2k+1`), then
    /// $$
    /// (2k+1)!! = (2k+1)(2k-1)(2k-3)\ldots(3)(1)
    /// $$
    /// Additionally we have that $0!!=1$.
    ///
    /// # Examples
    /// ## For `usize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0.factorial2(), 1);
    /// assert_eq!(3.factorial2(), 3);
    /// assert_eq!(5.factorial2(), 15);
    /// ```
    /// ## For `isize`
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(0_isize.factorial2(), 1);
    /// assert_eq!(3_isize.factorial2(), 3);
    /// assert_eq!(6_isize.factorial2(), 48);
    /// ```
    /// # Notes
    /// If `n < 0` returns `0`.
    fn factorial2(self) -> Self;

    /// Generalized `k`-factorial.
    ///
    /// Generalization of the [factorial] to steps of size `k`.
    /// $$
    /// n(!!\ldots!) = (n)(n-k)(n-2k)\ldots
    /// $$
    /// We always have that $0(!!\ldots!) = 1$ irregardless of `k`.
    ///
    /// In particular for any integer `n`, we have that
    /// ```
    /// # use sci_rs::special::Factorial;
    /// # let n = 1;
    ///
    /// assert_eq!(n.factorialk(1), n.factorial());
    /// assert_eq!(n.factorialk(1), n.factorial2());
    /// ```
    ///
    /// # Examples
    /// ```
    /// use sci_rs::special::Factorial;
    ///
    /// assert_eq!(5.factorialk(3), 10); // 5 * 2
    /// assert_eq!(10.factorialk(5), 50); // 10 * 5
    /// ```
    /// # Notes
    /// If `n < 0` returns `0`. `k` must be strictly greater than `0`.
    ///
    /// [factorial]: crate::special::Factorial::factorial
    fn factorialk(self, k: Self) -> Self;
}

macro_rules! factorial_primint_impl {
    ($($T: ty)*) => ($(
        impl Factorial for $T {
            #[inline(always)]
            fn factorial(self) -> Self {
                primint_factorial(self)
            }

            #[inline(always)]
            fn factorial2(self) -> Self {
                primint_factorial2(self)
            }

            #[inline(always)]
            fn factorialk(self, k: Self) -> Self {
                primint_factorialk(self, k)
            }
        }
    )*)
}

factorial_primint_impl! {u8 u16 u32 u64 u128 usize i8 i16 i32 i64 i128 isize}

const FACTORIAL_CACHE_LEN: usize = 17;
const FACTORIAL_CACHE: [u128; FACTORIAL_CACHE_LEN] = [
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
];

fn primint_factorial<Int>(n: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if n < Int::zero() {
        return Int::zero();
    }
    if n < Int::from_usize(FACTORIAL_CACHE_LEN).unwrap() {
        return Int::from_u128(FACTORIAL_CACHE[n.to_usize().unwrap()]).unwrap();
    }

    Int::from_u128(FACTORIAL_CACHE[FACTORIAL_CACHE_LEN - 1]).unwrap()
        * partial_product(Int::from_usize(FACTORIAL_CACHE_LEN).unwrap(), n, Int::one())
}

fn primint_factorial2<Int>(n: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    if n < Int::zero() {
        return Int::zero();
    }
    if n < Int::from_usize(FACTORIAL2_CACHE_LEN).unwrap() {
        return Int::from_u128(FACTORIAL2_CACHE[n.to_usize().unwrap()]).unwrap();
    }

    let two = Int::one() + Int::one();
    if (n & Int::one()).is_zero() && (FACTORIAL2_CACHE_LEN & 1 != 0) {
        Int::from_u128(FACTORIAL2_CACHE[FACTORIAL2_CACHE_LEN - 1]).unwrap()
            * partial_product(Int::from_usize(FACTORIAL2_CACHE_LEN + 1).unwrap(), n, two)
    } else {
        Int::from_u128(FACTORIAL2_CACHE[FACTORIAL2_CACHE_LEN - 2]).unwrap()
            * partial_product(Int::from_usize(FACTORIAL2_CACHE_LEN).unwrap(), n, two)
    }
}

const FACTORIAL2_CACHE_LEN: usize = 33;
const FACTORIAL2_CACHE: [u128; FACTORIAL2_CACHE_LEN] = [
    1,
    1,
    2,
    3,
    8,
    15,
    48,
    105,
    384,
    945,
    3840,
    10395,
    46080,
    135135,
    645120,
    2027025,
    10321920,
    34459425,
    185794560,
    654729075,
    3715891200,
    13749310575,
    81749606400,
    316234143225,
    1961990553600,
    7905853580625,
    51011754393600,
    213458046676875,
    1428329123020800,
    6190283353629375,
    42849873690624000,
    191898783962510625,
    1371195958099968000,
];

fn primint_factorialk<Int>(n: Int, k: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    assert!(k > Int::zero());

    if n < Int::zero() {
        return Int::zero();
    }

    if n.is_zero() {
        return Int::one();
    }

    let start = n - k * (n / k);
    if start.is_zero() {
        partial_product(k, n, k)
    } else {
        partial_product(start, n, k)
    }
}

/// Computes the product between `start` and `stop` stepping with `step`.
///
/// $$
/// n_0\times n_1\times\ldots \times n_k
/// $$
/// where $n_0 = $ `start`, $n_{i+1} = \text{step}\times n_i$, and $n_k$
/// is the largest value such that $n_k < $ `stop + 1`./
fn partial_product<Int>(mut start: Int, stop: Int, step: Int) -> Int
where
    Int: PrimInt + FromPrimitive,
{
    assert!(step > Int::zero());
    assert!(stop >= start);

    if step.is_one() {
        let start = start.to_usize().unwrap();
        let stop = stop.to_usize().unwrap();

        return (start..stop + 1)
            .map(|n| Int::from_usize(n).unwrap())
            .fold(Int::one(), |result, val| result * val);
    }

    let mut result = Int::one();

    while start <= stop {
        result = result * start;
        start = start + step;
    }
    result
}

#[cfg(test)]
mod test {
    use core::{any::type_name, fmt};
    use num_traits::ToPrimitive;

    use super::*;

    const ABSOLUTE_KNOWN_FACTORIAL_VALUES: [u128; 35] = [
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800,
        87178291200,
        1307674368000,
        20922789888000,
        355687428096000,
        6402373705728000,
        121645100408832000,
        2432902008176640000,
        51090942171709440000,
        1124000727777607680000,
        25852016738884976640000,
        620448401733239439360000,
        15511210043330985984000000,
        403291461126605635584000000,
        10888869450418352160768000000,
        304888344611713860501504000000,
        8841761993739701954543616000000,
        265252859812191058636308480000000,
        8222838654177922817725562880000000,
        263130836933693530167218012160000000,
        8683317618811886495518194401280000000,
        295232799039604140847618609643520000000,
    ];

    const ABSOLUTE_KNOWN_FACTORIAL2_VALUES: [u128; 57] = [
        1,
        1,
        2,
        3,
        8,
        15,
        48,
        105,
        384,
        945,
        3840,
        10395,
        46080,
        135135,
        645120,
        2027025,
        10321920,
        34459425,
        185794560,
        654729075,
        3715891200,
        13749310575,
        81749606400,
        316234143225,
        1961990553600,
        7905853580625,
        51011754393600,
        213458046676875,
        1428329123020800,
        6190283353629375,
        42849873690624000,
        191898783962510625,
        1371195958099968000,
        6332659870762850625,
        46620662575398912000,
        221643095476699771875,
        1678343852714360832000,
        8200794532637891559375,
        63777066403145711616000,
        319830986772877770815625,
        2551082656125828464640000,
        13113070457687988603440625,
        107145471557284795514880000,
        563862029680583509947946875,
        4714400748520531002654720000,
        25373791335626257947657609375,
        216862434431944426122117120000,
        1192568192774434123539907640625,
        10409396852733332453861621760000,
        58435841445947272053455474390625,
        520469842636666622693081088000000,
        2980227913743310874726229193921875,
        27064431817106664380040216576000000,
        157952079428395476360490147277859375,
        1461479318123759876522171695104000000,
        8687364368561751199826958100282265625,
        81842841814930553085241614925824000000,
    ];

    fn check_factorial<Int>(max: Int)
    where
        Int: Factorial + FromPrimitive + ToPrimitive + fmt::Debug + fmt::Display + PartialEq,
    {
        for index in 0..max.to_usize().unwrap() {
            let ref_value = Int::from_u128(ABSOLUTE_KNOWN_FACTORIAL_VALUES[index]).unwrap();
            assert_eq!(
                Int::factorial(Int::from_usize(index).unwrap()),
                ref_value,
                "{}! != {} (as {})",
                index,
                ref_value,
                type_name::<Int>()
            );
        }
    }

    fn check_factorial2<Int>(max: Int)
    where
        Int: Factorial + FromPrimitive + ToPrimitive + fmt::Debug + fmt::Display + PartialEq,
    {
        for index in 0..max.to_usize().unwrap() {
            let ref_value = Int::from_u128(ABSOLUTE_KNOWN_FACTORIAL2_VALUES[index]).unwrap();
            assert_eq!(
                Int::factorial2(Int::from_usize(index).unwrap()),
                ref_value,
                "{}!! != {} (as {})",
                index,
                ref_value,
                type_name::<Int>()
            );
        }
    }

    #[test]
    fn factorial() {
        check_factorial::<u8>(6);
        check_factorial::<u16>(9);
        check_factorial::<u32>(13);
        check_factorial::<u64>(21);
        check_factorial::<u128>(35);

        check_factorial::<i8>(6);
        check_factorial::<i16>(8);
        check_factorial::<i32>(13);
        check_factorial::<i64>(21);
        check_factorial::<i128>(34);
    }

    #[test]
    fn factorial2() {
        check_factorial2::<u8>(8);
        check_factorial2::<u16>(13);
        check_factorial2::<u32>(21);
        check_factorial2::<u64>(34);
        check_factorial2::<u128>(57);

        check_factorial2::<i8>(8);
        check_factorial2::<i16>(12);
        check_factorial2::<i32>(20);
        check_factorial2::<i64>(34);
        check_factorial2::<i128>(57);
    }

    #[test]
    fn factorialk() {
        const REF_VALUES: [[usize; 14]; 15] = [
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [6, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [24, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
            [120, 15, 10, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
            [720, 48, 18, 12, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
            [5040, 105, 28, 21, 14, 7, 7, 7, 7, 7, 7, 7, 7, 7],
            [40320, 384, 80, 32, 24, 16, 8, 8, 8, 8, 8, 8, 8, 8],
            [362880, 945, 162, 45, 36, 27, 18, 9, 9, 9, 9, 9, 9, 9],
            [
                3628800, 3840, 280, 120, 50, 40, 30, 20, 10, 10, 10, 10, 10, 10,
            ],
            [
                39916800, 10395, 880, 231, 66, 55, 44, 33, 22, 11, 11, 11, 11, 11,
            ],
            [
                479001600, 46080, 1944, 384, 168, 72, 60, 48, 36, 24, 12, 12, 12, 12,
            ],
            [
                6227020800, 135135, 3640, 585, 312, 91, 78, 65, 52, 39, 26, 13, 13, 13,
            ],
            [
                87178291200,
                645120,
                12320,
                1680,
                504,
                224,
                98,
                84,
                70,
                56,
                42,
                28,
                14,
                14,
            ],
        ];

        for (n, &elements) in REF_VALUES.iter().enumerate() {
            for (k, &value) in elements.iter().enumerate() {
                assert_eq!(n.factorialk(k + 1), value);
                assert_eq!(
                    n.to_i64().unwrap().factorialk((k + 1).to_i64().unwrap()),
                    value.to_i64().unwrap()
                );
            }
        }
    }
}
