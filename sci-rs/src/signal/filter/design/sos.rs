use itertools::Itertools;
use nalgebra::RealField;
use num_traits::Float;

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

    /// Filter delay value
    pub zi0: F,
    /// Filter delay value
    pub zi1: F,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(all(feature = "alloc", feature = "std"))]
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
