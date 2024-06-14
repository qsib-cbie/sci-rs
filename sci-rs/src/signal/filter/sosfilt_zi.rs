use core::{iter::Sum, ops::SubAssign};
use nalgebra::{ClosedAdd, ClosedMul, RealField, Scalar};
use num_traits::{Float, One, Zero};

use crate::signal::filter::lfilter_zi_dyn;

use super::design::Sos;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

///
/// Construct initial conditions for sosfilt for step response steady-state.
///
/// Compute an initial state `zi` for the `sosfilt` function that corresponds
/// to the steady state of the step response.
///
/// A typical use of this function is to set the initial state so that the
/// output of the filter starts at the same value as the first element of
/// the signal to be filtered.
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter_zi.html#scipy.signal.lfilter_zi>
///
///
pub fn sosfilt_zi_dyn<'a, F, I, S>(s: I)
where
    F: RealField + Copy + PartialEq + Scalar + Zero + One + ClosedMul + ClosedAdd + Sum + SubAssign,
    I: Iterator<Item = &'a mut Sos<F>>,
{
    let mut scale = F::one();
    for section in s {
        let zi = lfilter_zi_dyn(&section.b, &section.a);
        section.zi0 = scale * zi[0];
        section.zi1 = scale * zi[1];
        scale *= section.b.iter().cloned().sum::<F>() / section.a.iter().cloned().sum::<F>();
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use nalgebra::matrix;

    use super::*;

    #[test]
    fn scipy_example_dyn() {
        // sos = signal.butter(9, 0.125, output='sos')
        let mut sos: Vec<Sos<f64>> = Sos::from_scipy_dyn(
            5,
            [
                1.55854712e-07,
                3.11709423e-07,
                1.55854712e-07,
                1.00000000e+00,
                -6.68178638e-01,
                0.00000000e+00,
                1.00000000e+00,
                2.00000000e+00,
                1.00000000e+00,
                1.00000000e+00,
                -1.35904130e+00,
                4.71015698e-01,
                1.00000000e+00,
                2.00000000e+00,
                1.00000000e+00,
                1.00000000e+00,
                -1.42887946e+00,
                5.46607979e-01,
                1.00000000e+00,
                2.00000000e+00,
                1.00000000e+00,
                1.00000000e+00,
                -1.55098998e+00,
                6.78779458e-01,
                1.00000000e+00,
                1.00000000e+00,
                0.00000000e+00,
                1.00000000e+00,
                -1.73262236e+00,
                8.75376926e-01,
            ]
            .to_vec(),
        );
        assert_eq!(sos.len(), 5);

        // zi = signal.sosfilt_zi(sos)
        let expected_zi = matrix!(
                1.72292381e-06,  1.55854712e-07;
                6.52357932e-05, -2.97332383e-05;
                2.21320188e-03, -1.17932460e-03;
                6.90969677e-02, -4.61691178e-02;
                9.28622716e-01, -8.75376926e-01
        );

        // Compute zi inplace on Sos
        sosfilt_zi_dyn::<_, _, Sos<f64>>(sos.iter_mut());

        for (i, row) in expected_zi.row_iter().enumerate() {
            let section = sos[i];
            assert_relative_eq!(row[0], section.zi0, max_relative = 1e-6);
            assert_relative_eq!(row[1], section.zi1, max_relative = 1e-6);
        }
    }
}
