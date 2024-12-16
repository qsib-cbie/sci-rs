use nalgebra::RealField;
use ndarray::{Array, ArrayBase, Data, Dimension, RawData};

/// """
/// Return a periodic square-wave waveform.
///
/// The square wave has a period ``2*pi``, has value +1 from 0 to
/// ``2*pi*duty`` and -1 from ``2*pi*duty`` to ``2*pi``. `duty` must be in
/// the interval \[0,1\].
///
/// Note that this is not band-limited.  It produces an infinite number
/// of harmonics, which are aliased back and forth across the frequency
/// spectrum.
///
/// Parameters
/// ----------
/// t : array_like
///     The input time array.
/// duty : array_like, optional
///     Duty cycle.  Default is 0.5 (50% duty cycle).
///     If an array, causes wave shape to change over time, and must be the
///     same length as t.
///
/// Returns
/// -------
/// y : ndarray
///     Output array containing the square waveform.
///
/// Examples
/// --------
/// A 5 Hz waveform sampled at 500 Hz for 1 second:
///
/// >>> import numpy as np
/// >>> from scipy import signal
/// >>> import matplotlib.pyplot as plt
/// >>> t = np.linspace(0, 1, 500, endpoint=False)
/// >>> plt.plot(t, signal.square(2 * np.pi * 5 * t))
/// >>> plt.ylim(-2, 2)
///
/// A pulse-width modulated sine wave:
///
/// >>> plt.figure()
/// >>> sig = np.sin(2 * np.pi * t)
/// >>> pwm = signal.square(2 * np.pi * 30 * t, duty=(sig + 1)/2)
/// >>> plt.subplot(2, 1, 1)
/// >>> plt.plot(t, sig)
/// >>> plt.subplot(2, 1, 2)
/// >>> plt.plot(t, pwm)
/// >>> plt.ylim(-1.5, 1.5)
///
/// """
pub fn square<F, S, D>(t: &ArrayBase<S, D>, duty: F) -> Array<F, D>
where
    F: RealField,
    S: Data<Elem = F>,
    D: Dimension,
{
    assert!(F::zero() <= duty && duty <= F::one());
    let duty_threshold = F::two_pi() * duty;
    t.mapv(|t| {
        let x = t % F::two_pi();
        // Because % is the reminder and not the modulo operator, x can be negative.
        let x = if x < F::zero() { x + F::two_pi() } else { x };
        debug_assert!(F::zero() <= x && x <= F::two_pi());
        if x < duty_threshold {
            F::one()
        } else {
            -F::one()
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::{arr1, arr3};

    #[test]
    fn test_square_zero_duty() {
        let t = arr1(&[
            -4.821, -4.15, -3.394, -3.386, -2.966, -2.735, -2.464, -2.277, -2.094, -0.8963,
            0.03853, 1.432, 2.384, 2.522, 2.732, 3.125, 3.297, 3.517, 3.602, 4.908,
        ]);
        let expected = arr1(&[
            -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
            -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
        ]);
        let result = square(&t, 0.0);
        assert_vec_eq(result, expected);
    }

    #[test]
    fn test_square_one_duty() {
        let t = arr1(&[
            -3.521, -3.284, -3.257, -2.367, -1.965, -1.933, -0.5399, -0.4277, -0.3761, 0.3024,
            0.3624, 0.6161, 0.784, 1.42, 1.869, 1.934, 4.0, 4.298, 4.389, 4.807,
        ]);
        let expected = arr1(&[
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0,
        ]);
        let result = square(&t, 1.0);
        assert_vec_eq(result, expected);
    }

    #[test]
    fn test_square_duty_03_05_07() {
        let t = arr1(&[
            -4.991, -4.973, -3.988, -3.084, -2.562, -2.378, -1.618, -1.449, -0.8056, -0.6883,
            -0.5677, -0.5353, -0.1377, -0.1142, 0.8072, 0.821, 1.836, 2.722, 4.189, 4.384,
        ]);

        let expected_03 = arr1(&[
            1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0,
            1.0, 1.0, -1.0, -1.0, -1.0,
        ]);
        let expected_05 = arr1(&[
            1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0,
            1.0, 1.0, 1.0, -1.0, -1.0,
        ]);
        let expected_07 = arr1(&[
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0,
        ]);
        let result_03 = square(&t, 0.3);
        let result_05 = square(&t, 0.5);
        let result_07 = square(&t, 0.7);
        assert_vec_eq(result_03, expected_03);
        assert_vec_eq(result_05, expected_05);
        assert_vec_eq(result_07, expected_07);
    }

    #[test]
    fn test_square_3d() {
        let t = arr3(&[
            [
                [-4.452, -4.182, -3.663, -3.307, -2.995],
                [-2.482, -2.46, -1.929, -1.823, -1.44],
            ],
            [
                [-0.8743, 0.5359, 0.9073, 2.101, 2.161],
                [2.582, 2.977, 3.966, 4.298, 4.659],
            ],
        ]);
        let expected = arr3(&[
            [[1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, -1.0, -1.0, -1.0]],
            [[-1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, -1.0, -1.0]],
        ]);
        let result = square(&t, 0.67);
        assert_vec_eq(result, expected);
    }

    #[track_caller]
    fn assert_vec_eq<D: Dimension>(a: Array<f32, D>, b: Array<f32, D>) {
        assert_eq!(a.shape(), b.shape());
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
