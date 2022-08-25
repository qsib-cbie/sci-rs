use core::{f64::consts::PI, ops::Mul};

use nalgebra::{Complex, ComplexField, RealField};
use num_traits::Float;

use super::{FilterBandType, FilterOutput, FilterOutputType, FilterType, Sos, Zpk};

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
/// >>> plt.xlabel('Frequency [Hz]')
/// >>> plt.ylabel('Magnitude [dB]')
/// >>> plt.grid(True)
/// """
///
pub fn bilinear_zpk<F, const N: usize>(zpk: Zpk<F, N>, fs: F) -> Zpk<F, N>
where
    F: RealField + Float,
    [Sos<F>; N / 2 - 1]: Sized,
{
    todo!()
    // z = atleast_1d(z)
    // p = atleast_1d(p)

    // degree = _relative_degree(z, p)

    // fs2 = 2.0*fs

    // # Bilinear transform the poles and zeros
    // z_z = (fs2 + z) / (fs2 - z)
    // p_z = (fs2 + p) / (fs2 - p)

    // # Any zeros that were at infinity get moved to the Nyquist frequency
    // z_z = append(z_z, -ones(degree))

    // # Compensate for gain change
    // k_z = k * real(prod(fs2 - z) / prod(fs2 - p))

    // return z_z, p_z, k_z
}
