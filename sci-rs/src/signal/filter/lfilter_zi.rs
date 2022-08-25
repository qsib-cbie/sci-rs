use core::{iter::Sum, ops::SubAssign};
use heapless::Vec;
use nalgebra::{ClosedAdd, ClosedMul, DMatrix, RealField, SMatrix, Scalar};
use num_traits::{Float, One, Zero};

use crate::linalg::companion_dyn;

///
/// Construct initial conditions for lfilter for step response steady-state.
///
/// Compute an initial state `zi` for the `lfilter` function that corresponds
/// to the steady state of the step response.
///
/// A typical use of this function is to set the initial state so that the
/// output of the filter starts at the same value as the first element of
/// the signal to be filtered.
///
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.lfilter_zi.html#scipy.signal.lfilter_zi>
///
///
pub fn lfilter_zi_dyn<F, const M: usize>(b: &[F; M], a: &[F; M]) -> [F; M - 1]
where
    F: RealField + Copy + PartialEq + Scalar + Zero + One + ClosedMul + ClosedAdd + Sum + SubAssign,
{
    let ai0 = a
        .iter()
        .enumerate().find(|(_, ai)| **ai != F::zero())
        .expect("There must be at least one nonzero `a` coefficient.")
        .0;

    // Mormalize to a[0] == 1
    let mut a = a.iter().skip(ai0).cloned().collect::<Vec<_, M>>();
    let mut b = b.iter().cloned().collect::<Vec<_, M>>();
    let a0 = a[0];
    if a0 != F::one() {
        a = a.iter_mut().map(|xi| *xi / a0).collect();
        b = b.iter_mut().map(|xi| *xi / a0).collect();
    }

    // Pad with zeros to match length
    while a.len() < M {
        a.push(F::zero()).unwrap();
    }

    // Solve zi = A*zi + B
    let mut compa: DMatrix<_> = companion_dyn(a.iter(), M);
    compa.transpose_mut();
    let i_minus_a: SMatrix<F, { M - 1 }, { M - 1 }> = SMatrix::one() - compa;
    let mut zi: [F; M - 1] = [F::zero(); { M - 1 }];
    let b0 = b[0];
    let bsum = b
        .iter()
        .skip(1)
        .zip(a.iter().skip(1))
        .map(|(bi, ai)| *bi - *ai * b0)
        .sum::<F>();
    let z0 = bsum / i_minus_a.column(0).sum();
    zi[0] = z0;

    let mut asum = F::one();
    let mut csum = F::zero();
    let b0 = b[0];
    for k in 1..M - 1 {
        unsafe {
            let ak = *a.get_unchecked(k);
            asum += ak;
            csum += *b.get_unchecked(k) - ak * b0;
            *zi.get_unchecked_mut(k) = asum * z0 - csum;
        }
    }

    zi
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn scipy_example_dyn() {
        // b, a = butter(5, 0.25)
        let b = [
            0.00327922, 0.01639608, 0.03279216, 0.03279216, 0.01639608, 0.00327922,
        ];
        let a = [
            1.,
            -2.47441617,
            2.81100631,
            -1.70377224,
            0.54443269,
            -0.07231567,
        ];

        // zi = lfilter_zi(b, a)
        let expected_zi = [0.99672078, -1.49409147, 1.28412268, -0.45244173, 0.07559489];

        let zi = lfilter_zi_dyn(&b, &a);
        expected_zi.iter().zip(zi.iter()).for_each(|(e, a)| {
            assert_relative_eq!(e, a, max_relative = 1e-6);
        })
    }
}
