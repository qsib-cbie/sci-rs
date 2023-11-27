use itertools::Itertools;
use nalgebra::RealField;
use num_traits::Float;

use fixed::{
    traits::{FromFixed, ToFixed},
    types::extra::{U15, U28},
    FixedI16, FixedI32,
};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

///
/// Second Order Section Representation
/// of a digital filter transfer function
///
/// Example as Biquad Digital filter from
/// <https://en.wikipedia.org/wiki/Digital_biquad_filter>
///
/// H(z) =
///        (b0 + b1 * z_1 + b2 * z_2) /
///        (a0 + a1 * z_1 + a2 * z_2)
///
/// Sos { b: [b0, b1, b2], a: [a0, a1, a2] }
///
#[derive(Debug, Copy, Clone)]
pub struct Sos<F: RealField + Copy> {
    /// Transfer coefficients numerator
    pub b: [F; 3],
    /// Transfer coefficients denominator
    pub a: [F; 3],

    ///
    /// Filter Delay Values
    ///
    pub(crate) zi0: F,
    pub(crate) zi1: F,
}

impl<F: RealField + Copy> Default for Sos<F> {
    fn default() -> Self {
        Self {
            b: [F::zero(); 3],
            a: [F::zero(); 3],
            zi0: F::zero(),
            zi1: F::zero(),
        }
    }
}

impl<F: RealField + Copy> Sos<F> {
    /// Create a Second Order Section Biquad
    /// with the coefficients `b` and `a`. The
    /// coefficients need not be negated as they should match
    /// SciPy's output. Some other libraries like CMSIS DSP
    /// negate some coefficients in `a`.
    pub fn new(b: [F; 3], a: [F; 3]) -> Sos<F> {
        Sos::<F> {
            b,
            a,
            ..Default::default()
        }
    }

    /// Create a Second Order Section Biquad directly
    /// from a vector floats from scipy
    #[cfg(feature = "alloc")]
    pub fn from_scipy_dyn(order: usize, sos: Vec<F>) -> Vec<Sos<F>> {
        assert!(order * 6 == sos.len());

        sos.iter()
            .tuples::<(&F, &F, &F, &F, &F, &F)>()
            .take(order)
            .map(|ba| Sos::new([*ba.0, *ba.1, *ba.2], [*ba.3, *ba.4, *ba.5]))
            .collect()
    }
}

///
/// A fixed point Q4.24 number representing
/// a 32 bit signed integer with 24 fractional bits
///
/// [-1, 1) is represented by [-2^24, 2^24)
///
pub type Q28 = FixedI32<U28>;

///
/// A fixed point Q1.15 number representing
/// a 16 bit signed integer with 15 fractional bits
/// [-1, 1) is represented by [-2^15, 2^15)
///
pub type Q15 = FixedI16<U15>;

///
/// A second order section of Q8.24 fixed point coefficients and state
///
#[derive(Debug, Copy, Clone)]
pub struct SosQ28 {
    /// Transfer coefficients numerator
    pub b: [Q28; 3],
    /// Transfer coefficients denominator
    pub a: [Q28; 3],

    ///
    /// Filter Delay Values
    ///
    pub(crate) zi0: Q28,
    pub(crate) zi1: Q28,
}

impl SosQ28 {
    ///
    /// Convert a floating point Sos to a fixed point Sos
    /// with Q8.24 coefficients and state
    ///
    /// The coefficients must be scaled within the range of
    /// [-1, 1) to fit within the Q8.24 fixed point format
    ///
    #[cfg(feature = "alloc")]
    pub fn from_sos<F>(sos: &::alloc::vec::Vec<Sos<F>>) -> ::alloc::vec::Vec<SosQ28>
    where
        F: ToFixed + RealField + Copy,
    {
        let sections = sos
            .iter()
            .map(|sos| SosQ28 {
                b: [
                    Q28::from_num((sos.b[0])),
                    Q28::from_num((sos.b[1])),
                    Q28::from_num((sos.b[2])),
                ],
                a: [
                    Q28::from_num((sos.a[0])),
                    Q28::from_num((sos.a[1])),
                    Q28::from_num((sos.a[2])),
                ],
                zi0: Q28::ZERO,
                zi1: Q28::ZERO,
            })
            .collect();
        sections
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "alloc")]
    #[test]
    fn can_use_external_filter_design() {
        let _design_filter = r#"
import scipy.signal as sg
buttersos = sg.butter(4, [0.5, 100], btype='bandpass', output='sos', fs=1666)
print(buttersos)
        "#;

        let butterworth_filter_sos = [
            0.00474269,
            0.00948539,
            0.00474269,
            1.,
            -1.05531479,
            0.29986557,
            1.,
            2.,
            1.,
            1.,
            -1.32397785,
            0.6355536,
            1.,
            -2.,
            1.,
            1.,
            -1.99416225,
            0.99417226,
            1.,
            -2.,
            1.,
            1.,
            -1.99760501,
            0.9976149,
        ];
        let sos: Vec<Sos<f64>> = Sos::from_scipy_dyn(4, butterworth_filter_sos.to_vec());
        assert_eq!(sos.len(), 4);
        println!("{:?}", sos);
    }
}
