use itertools::Itertools;
use nalgebra::RealField;
use num_traits::Float;

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
    ///
    /// Transfer coefficients
    ///
    pub b: [F; 3],
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
    // Normalizes to Direct Form 1 where a[0] is 1
    pub fn new(mut b: [F; 3], mut a: [F; 3]) -> Sos<F> {
        b[0] = b[0] / a[0];
        b[1] = b[1] / a[0];
        b[2] = b[2] / a[0];
        a[0] = F::one();
        a[1] = a[1] / a[0];
        a[2] = a[2] / a[0];
        Sos::<F> {
            b,
            a,
            ..Default::default()
        }
    }

    pub fn from_scipy<const M: usize, const N: usize>(sos: [F; M]) -> [Sos<F>; N] {
        // TODO: Replace with stable const expressions
        assert!(N * 6 == M);

        let mut rslt = [Default::default(); N];
        sos.iter()
            .tuples::<(&F, &F, &F, &F, &F, &F)>()
            .take(N)
            .zip(rslt.iter_mut())
            .for_each(|(ba, sos)| {
                *sos = Sos::new([*ba.0, *ba.1, *ba.2], [*ba.3, *ba.4, *ba.5]);
            });
        rslt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let sos: [Sos<f64>; 4] = Sos::from_scipy(butterworth_filter_sos);
        println!("{:?}", sos);
    }
}
