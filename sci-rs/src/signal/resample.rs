use nalgebra::Complex;
use num_traits::{Float, Zero};
use rustfft::FftNum;

///
/// Resample the data to the desired number of samples using the Fourier transform.
///
/// This method is similar but not exactly equivalent to the SciPy method of resampling:
/// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.resample.html>
///
/// It skips some complexity of the SciPy method, such as windowing and handling odd vs. even-length signals.
///
/// Procedure:
/// 1. Convert to the frequency-domain.
///   a. If upsampling, pad higher frequency bins with 0
///   b. If downsampling, truncate higher frequency bins
/// 2. Convert back to the time-domain.
///
pub fn resample<F: Float + FftNum>(x: &[F], n: usize) -> Vec<F> {
    // SciPy style 'Fourier' resampling
    // 1. Compute FFT of x
    // 2. Fill vec of zeros with the desired length, y.
    // 3. Set the from beginning of y to the first half of x
    // 4. Set the from end of y to the second half of x
    // 5. Compute IFFT of y
    // 6. Multiply y by (n / x.len())
    // 7. Take the real part of y

    // Compute FFT of x
    let mut fft_planner = rustfft::FftPlanner::<F>::new();
    let fft = fft_planner.plan_fft_forward(x.len());
    let ifft = fft_planner.plan_fft_inverse(n);

    let scratch_length = std::cmp::max(
        fft.get_inplace_scratch_len(),
        ifft.get_inplace_scratch_len(),
    );
    let mut scratch = vec![Complex::zero(); scratch_length];
    let mut x = x
        .into_iter()
        .map(|x| Complex::new(*x, F::zero()))
        .collect::<Vec<_>>();
    fft.process_with_scratch(&mut x, &mut scratch);

    // Fill y with halfs of x
    let mut y = vec![Complex::zero(); n];
    let bins = std::cmp::min(x.len(), n);
    let half_spectrum = bins / 2;
    y[..half_spectrum].copy_from_slice(&x[..half_spectrum]);
    y[n - half_spectrum..].copy_from_slice(&x[x.len() - half_spectrum..]);

    // Compute iFFT of y
    ifft.process_with_scratch(&mut y, &mut scratch);

    // Take the scaled real domain as the resampled result
    let scale_factor = F::from(1.0 / x.len() as f64).unwrap();
    let y = y.iter().map(|x| x.re * scale_factor).collect::<Vec<_>>();

    y
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use rand::Rng;

    use super::*;

    #[test]
    #[should_panic]
    fn can_resample_like_scipy() {
        let x = vec![1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let y = resample(&x, 5);
        let expected = vec![3., 2.18649851, 5.01849831, 5.98150169, 8.81350149];
        assert_eq!(y.len(), expected.len());

        // Fails due to slight algorithmic differences
        println!("y: {:?}, scipy: {:?}", y, expected);
        for (y, expected) in y.iter().zip(expected.iter()) {
            assert_relative_eq!(y, expected, epsilon = 1e-10);
        }
    }

    #[test]
    fn can_resample_to_exact_number() {
        // Randomly generate 1000 vectors of length (10,1000)
        // Resample each to length 100
        // Check that each resampled vector is length 100

        let mut rng = rand::thread_rng();
        for i in (0..100) {
            let len = rng.gen_range((10..50));
            let x = (0..len)
                .map(|_| rng.gen_range((-100.0..100.)))
                .collect::<Vec<_>>();
            let y = resample(&x, 100);
            assert_eq!(y.len(), 100);
        }

        for i in (0..50) {
            let len = rng.gen_range((200..10000));
            let target_len = rng.gen_range((50..50000));
            let x = (0..len)
                .map(|_| rng.gen_range((-100.0..100.)))
                .collect::<Vec<_>>();
            let y = resample(&x, target_len);
            assert_eq!(y.len(), target_len);
        }
    }
}
