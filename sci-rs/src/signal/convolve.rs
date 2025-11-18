use ndarray::{
    Array, ArrayBase, ArrayView, Data, Dim, IntoDimension, Ix, RemoveAxis, SliceArg, SliceInfo,
    SliceInfoElem,
};
use ndarray_conv::{ConvFFTExt, ConvMode};
use num_traits::NumAssign;
use rustfft::FftNum;

/// Convolution mode determines behavior near edges and output size
pub enum ConvolveMode {
    /// Full convolution, output size is `in1.len() + in2.len() - 1`
    Full,
    /// Valid convolution, output size is `max(in1.len(), in2.len()) - min(in1.len(), in2.len()) + 1`
    Valid,
    /// Same convolution, output size is `in1.len()`
    Same,
}

/// Convolve two N-dimensional arrays using the fourier method.
///
/// According to Python docs, this is generally much faster than direct convolution
/// for large arrays (n > ~500), but can be slower when only a few output values are needed.
///
/// # Arguments
/// - `in1`: First input signal by reference. Can be `[std::vec::Vec]` or `[ndarray::Array]`.
/// - `in2`: Second input signal by reference. (Same type and dimensions as `in1`.)
/// - `mode`: [ConvolveMode]
///
/// # Returns
/// An `[Array]` containing the discrete linear convolution of `in1` with `in2`.
/// For [ConvolveMode::Full] mode, the output length will be `in1.shape() "+" in2.shape() "-" 1`.
/// For [ConvolveMode::Valid] mode, the output length will be `max(in1.shape(), + in2.shape())`.
/// For [ConvolveMode::Same] mode, the output length will be `in1.shape()`.
pub fn fftconvolve<'a, T, S, const N: usize>(
    in1: ArrayBase<S, Dim<[Ix; N]>>,
    in2: ArrayBase<S, Dim<[Ix; N]>>,
    mode: ConvolveMode,
) -> Array<T, Dim<[Ix; N]>>
where
    T: NumAssign + FftNum,
    S: Data<Elem = T> + 'a,
    [Ix; N]: IntoDimension<Dim = Dim<[Ix; N]>>,
    Dim<[Ix; N]>: RemoveAxis,
    SliceInfo<[SliceInfoElem; N], Dim<[Ix; N]>, Dim<[Ix; N]>>:
        SliceArg<Dim<[Ix; N]>, OutDim = Dim<[Ix; N]>>,
{
    match mode {
        ConvolveMode::Full => {
            todo!()
        }
        ConvolveMode::Valid => {
            todo!()
        }
        ConvolveMode::Same => {
            in1.conv_fft(&in2, ConvMode::Same, ndarray_conv::PaddingMode::Zeros)
                .unwrap() // TODO: Result type from core
        }
    }
}

/// Compute the convolution of two signals using FFT.
///
/// # Arguments
/// - `in1`: First input signal by reference. Can be `[ndarray::Array]`.
/// - `in2`: Second input signal by reference. (Same type and dimensions as `in1`.)
/// - `mode`: [ConvolveMode]
///
/// # Returns
/// An `[Array]` containing the discrete linear convolution of `in1` with `in2`.
/// For [ConvolveMode::Full] mode, the output length will be `in1.shape() "+" in2.shape() "-" 1`.
/// For [ConvolveMode::Valid] mode, the output length will be `max(in1.shape(), + in2.shape())`.
/// For [ConvolveMode::Same] mode, the output length will be `in1.shape()`.
///
/// # Note
/// Automatic choice between convolution through direct summation or via FFT has yet to be done
#[inline]
pub fn convolve<'a, T, S, const N: usize>(
    in1: ArrayBase<S, Dim<[Ix; N]>>,
    in2: ArrayBase<S, Dim<[Ix; N]>>,
    mode: ConvolveMode,
) -> Array<T, Dim<[Ix; N]>>
where
    T: NumAssign + FftNum,
    S: Data<Elem = T> + 'a,
    [Ix; N]: IntoDimension<Dim = Dim<[Ix; N]>>,
    Dim<[Ix; N]>: RemoveAxis,
    SliceInfo<[SliceInfoElem; N], Dim<[Ix; N]>, Dim<[Ix; N]>>:
        SliceArg<Dim<[Ix; N]>, OutDim = Dim<[Ix; N]>>,
{
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
/// An array containing the cross-correlation of `in1` with `in2`.
/// With Full mode, the output length will be `in1.len() + in2.len() - 1`.
pub fn correlate<'a, T, S, const N: usize>(
    in1: ArrayBase<S, Dim<[Ix; N]>>,
    in2: ArrayBase<S, Dim<[Ix; N]>>,
    mode: ConvolveMode,
) -> Array<T, Dim<[Ix; N]>>
where
    T: NumAssign + FftNum,
    S: Data<Elem = T> + 'a,
    [Ix; N]: IntoDimension<Dim = Dim<[Ix; N]>>,
    Dim<[Ix; N]>: RemoveAxis,
    SliceInfo<[SliceInfoElem; N], Dim<[Ix; N]>, Dim<[Ix; N]>>:
        SliceArg<Dim<[Ix; N]>, OutDim = Dim<[Ix; N]>>,
{
    in1.conv_fft(
        &in2.t(),
        match mode {
            ConvolveMode::Full => ConvMode::Full,
            ConvolveMode::Valid => ConvMode::Valid,
            ConvolveMode::Same => ConvMode::Same,
        },
        ndarray_conv::PaddingMode::Zeros,
    )
    .unwrap() // TODO: Result type from core
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::{array, Array1, ArrayView1};

    #[test]
    fn test_convolve_full() {
        let in1 = array![1.0, 2.0, 3.0];
        let in2 = array![4.0, 5.0, 6.0];
        let result: Array1<_> = convolve(in1, in2, ConvolveMode::Full);
        let expected: Array1<_> = vec![4.0, 13.0, 28.0, 27.0, 18.0].into();

        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_correlate_full() {
        let in1 = array![1.0, 2.0, 3.0];
        let in2 = array![4.0, 5.0, 6.0];

        {
            let in1_ref: ArrayView1<_> = (&in1).into();
            let in2_ref: ArrayView1<_> = (&in2).into();
            let result: Array<f64, Dim<[Ix; 1]>> = correlate(in1_ref, in2_ref, ConvolveMode::Full);
            let expected: Array<f64, Dim<[Ix; 1]>> = vec![6.0, 17.0, 32.0, 23.0, 12.0].into();
            for (a, b) in result.iter().zip(expected.iter()) {
                assert_relative_eq!(a, b, epsilon = 1e-10);
            }
        }

        let result: Array1<_> = correlate(in1, in2, ConvolveMode::Full);
        let expected: Array1<_> = array![6.0, 17.0, 32.0, 23.0, 12.0];
        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_convolve_valid() {
        let in1 = array![1.0, 2.0, 5.0, 7.0];
        let in2 = array![1.4, 2.2];
        let result: Array1<_> = convolve(in1, in2, ConvolveMode::Valid);
        let expected: Array1<_> = array![5.0, 11.4, 20.8];
        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_convolve_same() {
        let in1 = array![1.0, 2.0, 3.0, 4.0];
        let in2 = array![1.0, 2.0, 1.5];
        let result: Array1<_> = convolve(in1, in2, ConvolveMode::Same);
        let expected: Array1<_> = array![4.0, 8.5, 13.0, 12.5];
        for (a, b) in result.iter().zip(expected.iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-10);
        }
    }

    #[test]
    #[cfg(feature = "plot")]
    fn test_scipy_example() {
        use rand::distributions::{Distribution, Standard};
        use rand::thread_rng;

        // Generate 1000 random samples from standard normal distribution
        let mut rng = thread_rng();
        let sig: Vec<f64> = Standard.sample_iter(&mut rng).take(1000).collect();
        let sig_ref: ArrayView1<_> = (&sig).into();

        // Compute autocorrelation using correlate directly
        let autocorr = correlate(sig_ref, sig_ref, ConvolveMode::Full);

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
