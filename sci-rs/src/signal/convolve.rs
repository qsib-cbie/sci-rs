use nalgebra::Complex;
use num_traits::{Float, FromPrimitive, Signed, Zero};
use rustfft::{FftNum, FftPlanner};

/// Convolution mode determines behavior near edges and output size
pub enum ConvolveMode {
    /// Full convolution, output size is `in1.len() + in2.len() - 1`
    Full,
    // Not yet implemented
    // Valid,
    // Same,
}

/// Performs FFT-based convolution on two slices of floating point values.
///
/// This is generally much faster than direct convolution for large arrays (n > ~500),
/// but can be slower when only a few output values are needed.
///
/// # Arguments
/// - `in1`: First input signal
/// - `in2`: Second input signal
/// - `mode`: Convolution mode (currently only Full is supported)
///
/// # Returns
/// A Vec containing the discrete linear convolution of `in1` with `in2`.
/// For Full mode, the output length will be `in1.len() + in2.len() - 1`.
pub fn fftconvolve<F: Float + FftNum>(in1: &[F], in2: &[F], mode: ConvolveMode) -> Vec<F> {
    match mode {
        ConvolveMode::Full => {
            // Determine the size of the FFT (next power of 2 for zero-padding)
            let n = in1.len() + in2.len() - 1;
            let fft_size = n.next_power_of_two();

            // Prepare input buffers as Complex<F> with zero-padding to fft_size
            let mut padded_in1 = vec![Complex::zero(); fft_size];
            let mut padded_in2 = vec![Complex::zero(); fft_size];

            // Copy input data into zero-padded buffers
            padded_in1.iter_mut().zip(in1.iter()).for_each(|(p, &v)| {
                *p = Complex::new(v, F::zero());
            });
            padded_in2.iter_mut().zip(in2.iter()).for_each(|(p, &v)| {
                *p = Complex::new(v, F::zero());
            });

            // Perform the FFT
            let mut planner = FftPlanner::new();
            let fft = planner.plan_fft_forward(fft_size);
            fft.process(&mut padded_in1);
            fft.process(&mut padded_in2);

            // Multiply element-wise in the frequency domain
            let mut result_freq: Vec<Complex<F>> = padded_in1
                .iter()
                .zip(&padded_in2)
                .map(|(a, b)| a * b)
                .collect();

            // Perform the inverse FFT
            let ifft = planner.plan_fft_inverse(fft_size);
            ifft.process(&mut result_freq);

            // Take only the real part, normalize, and truncate to the original output size (n)
            let fft_size = F::from(fft_size).unwrap();
            result_freq
                .iter()
                .take(n)
                .map(|x| x.re / fft_size)
                .collect()
        }
    }
}

/// Compute the convolution of two signals using FFT.
///
/// # Arguments
/// * `in1` - First input array
/// * `in2` - Second input array
///
/// # Returns
/// A Vec containing the convolution of `in1` with `in2`.
/// With Full mode, the output length will be `in1.len() + in2.len() - 1`.
pub fn convolve<F: Float + FftNum>(in1: &[F], in2: &[F], mode: ConvolveMode) -> Vec<F> {
    fftconvolve(in1, in2, mode)
}

/// Compute the cross-correlation of two signals using FFT.
///
/// Cross-correlation is similar to convolution but with flipping one of the signals.
/// This function uses FFT to compute the correlation efficiently.
///
/// # Arguments
/// * `in1` - First input array
/// * `in2` - Second input array
///
/// # Returns
/// A Vec containing the cross-correlation of `in1` with `in2`.
/// With Full mode, the output length will be `in1.len() + in2.len() - 1`.
pub fn correlate<F: Float + FftNum>(in1: &[F], in2: &[F], mode: ConvolveMode) -> Vec<F> {
    // For correlation, we need to reverse in2
    let mut in2_rev = in2.to_vec();
    in2_rev.reverse();
    fftconvolve(in1, &in2_rev, mode)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_convolve() {
        let in1 = vec![1.0, 2.0, 3.0];
        let in2 = vec![4.0, 5.0, 6.0];
        let result = convolve(&in1, &in2, ConvolveMode::Full);
        let expected = vec![4.0, 13.0, 28.0, 27.0, 18.0];

        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_correlate() {
        let in1 = vec![1.0, 2.0, 3.0];
        let in2 = vec![4.0, 5.0, 6.0];
        let result = correlate(&in1, &in2, ConvolveMode::Full);
        let expected = vec![6.0, 17.0, 32.0, 23.0, 12.0];
        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_scipy_example() {
        use rand::distributions::{Distribution, Standard};
        use rand::thread_rng;

        // Generate 1000 random samples from standard normal distribution
        let mut rng = thread_rng();
        let sig: Vec<f64> = Standard.sample_iter(&mut rng).take(1000).collect();

        // Compute autocorrelation using correlate directly
        let autocorr = correlate(&sig, &sig, ConvolveMode::Full);

        // Basic sanity checks
        assert_eq!(autocorr.len(), 1999); // Full convolution length should be 2N-1
        assert!(autocorr.iter().all(|&x| !x.is_nan())); // No NaN values

        // Maximum correlation should be near the middle since it's autocorrelation
        let max_idx = autocorr
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert!((max_idx as i32 - 999).abs() <= 1); // Should be near index 999

        let sig: Vec<f32> = sig.iter().map(|x| *x as f32).collect();
        let autocorr: Vec<f32> = autocorr.iter().map(|x| *x as f32).collect();
        crate::plot::python_plot(vec![&sig, &autocorr]);
    }
}
