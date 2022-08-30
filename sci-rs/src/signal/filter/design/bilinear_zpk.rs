use core::{f64::consts::PI, ops::Mul};

use heapless::Vec;
use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{
    relative_degree::relative_degree, FilterBandType, FilterOutputType, FilterType, Sos,
    ZpkFormatFilter,
};

/// """
/// Return a digital IIR filter from an analog one using a bilinear transform.
///
/// Transform a set of poles and zeros from the analog s-plane to the digital
/// z-plane using Tustin's method, which substitutes ``(z-1) / (z+1)`` for
/// ``s``, maintaining the shape of the frequency response.
///
/// Parameters
/// ----------
/// z : array_like
///     Zeros of the analog filter transfer function.
/// p : array_like
///     Poles of the analog filter transfer function.
/// k : float
///     System gain of the analog filter transfer function.
/// fs : float
///     Sample rate, as ordinary frequency (e.g., hertz). No prewarping is
///     done in this function.
///
/// Returns
/// -------
/// z : ndarray
///     Zeros of the transformed digital filter transfer function.
/// p : ndarray
///     Poles of the transformed digital filter transfer function.
/// k : float
///     System gain of the transformed digital filter.
///
/// See Also
/// --------
/// lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk
/// bilinear
///
/// Notes
/// -----
/// .. versionadded:: 1.1.0
///
/// Examples
/// --------
/// >>> from scipy import signal
/// >>> import matplotlib.pyplot as plt
///
/// >>> fs = 100
/// >>> bf = 2 * np.pi * np.array([7, 13])
/// >>> filts = signal.lti(*signal.butter(4, bf, btype='bandpass', analog=True,
/// ...                                   output='zpk'))
/// >>> filtz = signal.lti(*signal.bilinear_zpk(filts.zeros, filts.poles,
/// ...                                         filts.gain, fs))
/// >>> wz, hz = signal.freqz_zpk(filtz.zeros, filtz.poles, filtz.gain)
/// >>> ws, hs = signal.freqs_zpk(filts.zeros, filts.poles, filts.gain,
/// ...                           worN=fs*wz)
/// >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hz).clip(1e-15)),
/// ...              label=r'$|H_z(e^{j \omega})|$')
/// >>> plt.semilogx(wz*fs/(2*np.pi), 20*np.log10(np.abs(hs).clip(1e-15)),
/// ...              label=r'$|H(j \omega)|$')
/// >>> plt.legend()
/// >>> plt.xlabel('Frequency (Hz)')
/// >>> plt.ylabel('Magnitude (dB)')
/// >>> plt.grid(True)
/// """
///
pub fn bilinear_zpk<F, const N: usize>(zpk: ZpkFormatFilter<F, N>, fs: F) -> ZpkFormatFilter<F, N>
where
    F: RealField + Float,
{
    let degree = relative_degree(&zpk.z, &zpk.p);

    // Bilinear transform the poles and zeros
    let fs2 = Complex::new(F::from(2.).unwrap() * fs, F::zero());

    // # Any zeros that were at infinity get moved to the Nyquist frequency
    let z_z: Vec<Complex<F>, N> = zpk
        .z
        .iter()
        .map(|zi| (fs2 + zi) / (fs2 - zi))
        .chain((0..degree).map(|_| Complex::new(-F::one(), F::zero())))
        .collect::<Vec<_, N>>();
    let p_z: Vec<Complex<F>, N> = zpk
        .p
        .iter()
        .map(|pi| (fs2 + pi) / (fs2 - pi))
        .collect::<Vec<_, N>>();

    // Compensate for gain change
    let num = zpk
        .z
        .iter()
        .map(|zi| fs2 - *zi)
        .fold(Complex::new(F::one(), F::zero()), |acc, zi| acc * zi);
    let denom = zpk
        .p
        .iter()
        .map(|pi| fs2 - *pi)
        .fold(Complex::new(F::one(), F::zero()), |acc, pi| acc * pi);
    let k_z = zpk.k * (num / denom).real();

    ZpkFormatFilter::new(z_z, p_z, k_z)
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;
    use num_traits::Zero;

    use super::*;

    #[test]
    fn matches_scipy_iirfilter_butter() {
        // butter(4, [10, 50], btype='bandpass', output='sos', fs=1666)
        let fs = 2.;
        let zpk: ZpkFormatFilter<_, 8> = ZpkFormatFilter::new(
            Vec::from_slice(&[Complex::zero(); 4]).unwrap(),
            Vec::from_slice(&[
                Complex::new(-0.02022036, -0.07498294),
                Complex::new(-0.07648538, -0.06990013),
                Complex::new(-0.07648538, 0.06990013),
                Complex::new(-0.02022036, 0.07498294),
                Complex::new(-0.0956662, 0.35475786),
                Complex::new(-0.20328954, 0.1857867),
                Complex::new(-0.20328954, -0.1857867),
                Complex::new(-0.0956662, -0.35475786),
            ])
            .unwrap(),
            0.008409569194994788,
        );

        let expected_zpk: ZpkFormatFilter<_, 8> = ZpkFormatFilter::new(
            Vec::from_slice(&[
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
                Complex::new(-1., 0.),
            ])
            .unwrap(),
            Vec::from_slice(&[
                Complex::new(0.98924866, -0.03710237),
                Complex::new(0.96189799, -0.03364097),
                Complex::new(0.96189799, 0.03364097),
                Complex::new(0.98924866, 0.03710237),
                Complex::new(0.93873849, 0.16792939),
                Complex::new(0.89956011, 0.08396115),
                Complex::new(0.89956011, -0.08396115),
                Complex::new(0.93873849, -0.16792939),
            ])
            .unwrap(),
            2.6775767382597835e-05,
        );

        let actual_zpk: ZpkFormatFilter<f64, 8> = bilinear_zpk(zpk, fs);
        assert_eq!(actual_zpk.z.len(), expected_zpk.z.len());
        for (a, e) in actual_zpk.z.iter().zip(expected_zpk.z.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_eq!(actual_zpk.p.len(), expected_zpk.p.len());
        for (a, e) in actual_zpk.p.iter().zip(expected_zpk.p.iter()) {
            assert_relative_eq!(a.re, e.re, max_relative = 1e-6);
            assert_relative_eq!(a.im, e.im, max_relative = 1e-6);
        }

        assert_relative_eq!(actual_zpk.k, expected_zpk.k, max_relative = 1e-8);
    }
}
