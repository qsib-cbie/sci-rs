//! This file provides the modified Bessel function of order zero of the argument.
//! It has been written with reference to Numerical Recipes in C++ by William H. Press,
//! Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery.
//! The function as provisioned here are not *yet* intended to be exposed out of the crate.
//! See: https://github.com/scipy/scipy/pull/21297, scipy/special/xsf/cephes/i0.h for a different
//! implementation.

use num_traits::real::Real;

/// All [functions located in the `Faster versions of common Bessel
/// functions.`](<https://docs.scipy.org/doc/scipy/reference/special.html#faster-versions-of-common-bessel-functions>)
/// Note strongly that Bessel has tremendous inacurracy for f32 compared to f64.
pub(crate) trait Bessel {
    /// Modified Bessel function of order 0.
    ///
    /// ## Notes
    /// * The range is partitioned into the two intervals [0, 8] and (8, infinity).
    /// * [Scipy has this as a
    ///   ufunc](<https://docs.scipy.org/doc/scipy/reference/special.html#special-functions-scipy-special>),
    ///   as a supposed wrapper over the Cephes routine. We try to define it over reasonable types in
    ///   the impl.
    fn i0(&self) -> Self;
}

impl Bessel for f64 {
    fn i0(&self) -> Self {
        let x = *self;
        let ax = x.abs();
        match ax < 15.0 {
            true => unsafe { poly(&I0P, x * x) / poly(&I0Q, 225. - (x * x)) },
            false => unsafe {
                ax.exp() * poly(&I0PP, 1.0 - 15.0 / ax) / (poly(&I0QQ, 1.0 - 15.0 / ax) * ax.sqrt())
            },
        }
    }
}

impl Bessel for f32 {
    // Known to yield wrong result.
    fn i0(&self) -> Self {
        let x = *self;
        let ax = x.abs();
        match ax < 15.0 {
            true => unsafe { poly(&I0P_F32, x * x) / poly(&I0Q_F32, 225. - (x * x)) },
            false => unsafe {
                ax.exp() * poly(&I0PP_F32, 1.0 - 15.0 / ax)
                    / (poly(&I0QQ_F32, 1.0 - 15.0 / ax) * ax.sqrt())
            },
        }
    }
}

/// Evaluates a polynomial with coeffecients evaluated at `x`.
unsafe fn poly<T>(cof: &[T], x: T) -> T
where
    T: Real + core::ops::Mul,
{
    cof.iter()
        .take(cof.len() - 1)
        .rev()
        // SAFETY: poly is called only by const arrays
        .fold(unsafe { *cof.last().unwrap_unchecked() }, |acc, &e| {
            acc * x + e
        })
}

const I0P: [f64; 14] = [
    9.999999999999997e-1,
    2.466405579426905e-1,
    1.478980363444585e-2,
    3.826993559940360e-4,
    5.395676869878828e-6,
    4.700912200921704e-8,
    2.733894920915608e-10,
    1.115830108455192e-12,
    3.301093025084127e-15,
    7.209167098020555e-18,
    1.166898488777214e-20,
    1.378948246502109e-23,
    1.124884061857506e-26,
    5.498556929587117e-30,
];

const I0P_F32: [f32; 14] = [
    9.999999999999997e-1,
    2.466405579426905e-1,
    1.478980363444585e-2,
    3.826993559940360e-4,
    5.395676869878828e-6,
    4.700912200921704e-8,
    2.733894920915608e-10,
    1.115830108455192e-12,
    3.301093025084127e-15,
    7.209167098020555e-18,
    1.166898488777214e-20,
    1.378948246502109e-23,
    1.124884061857506e-26,
    5.498556929587117e-30,
];

const I0Q: [f64; 5] = [
    4.463598170691436e-1,
    1.702205745042606e-3,
    2.792125684538934e-6,
    2.369902034785866e-9,
    8.965900179621208e-13,
];

const I0Q_F32: [f32; 5] = [
    4.463598170691436e-1,
    1.702205745042606e-3,
    2.792125684538934e-6,
    2.369902034785866e-9,
    8.965900179621208e-13,
];

const I0PP: [f64; 5] = [
    1.192273748120670e-1,
    1.947452015979746e-1,
    7.629241821600588e-2,
    8.474903580801549e-3,
    2.023821945835647e-4,
];

const I0PP_F32: [f32; 5] = [
    1.192273748120670e-1,
    1.947452015979746e-1,
    7.629241821600588e-2,
    8.474903580801549e-3,
    2.023821945835647e-4,
];

const I0QQ: [f64; 6] = [
    2.962898424533095e-1,
    4.866115913196384e-1,
    1.938352806477617e-1,
    2.261671093400046e-2,
    6.450448095075585e-4,
    1.529835782400450e-6,
];

const I0QQ_F32: [f32; 6] = [
    2.962898424533095e-1,
    4.866115913196384e-1,
    1.938352806477617e-1,
    2.261671093400046e-2,
    6.450448095075585e-4,
    1.529835782400450e-6,
];

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn i0_f64() {
        let result: f64 = (1.).i0();
        let exp = 1.2660658777520082;
        assert_relative_eq!(result, exp, epsilon = 1e-6);
        let result: f64 = 0.213.i0();
        let exp = 1.0113744522192416;
        assert_relative_eq!(result, exp, epsilon = 1e-6);
        let result: f64 = 30.546.i0();
        let exp = 1337209608661.4026;
        assert_relative_eq!(result, exp, epsilon = 1e-6);
    }

    #[test]
    fn i0_f32() {
        let result: f32 = (1.).i0();
        let exp = 1.2660658777520082;
        assert_relative_eq!(result, exp, epsilon = 1e-6);
        let result: f32 = 0.213.i0();
        let exp = 1.0113744522192416;
        assert_relative_eq!(result, exp, epsilon = 1e-6);
    }
}
