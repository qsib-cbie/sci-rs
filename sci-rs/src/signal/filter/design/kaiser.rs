use nalgebra::RealField;
use num_traits::{real::Real, MulAdd, Pow, ToPrimitive};

/// Compute the Kaiser parameter `beta`, given the attenuation `a`.
///
/// # Parameters
/// * `a`: float  
///     The desired attenuation in the stopband and maximum ripple in the passband, in dB.  This
///     should be a *positive* number.
///
/// # Returns
/// * `beta`: float  
///     The `beta` parameter to be used in the formula for a Kaiser window.
///
/// # References
/// Oppenheim, Schafer, "Discrete-Time Signal Processing", p.475-476.
///
/// Examples
/// --------
/// Suppose we want to design a lowpass filter, with 65 dB attenuation
/// in the stop band.  The Kaiser window parameter to be used in the
/// window method is computed by ``kaiser_beta(65.)``:
///
/// ```
/// use sci_rs::signal::filter::design::kaiser_beta;
/// assert_eq!(6.20426, kaiser_beta(65.));
/// ```
/// # See Also
/// [kaiser_atten], [kaiserord]
pub fn kaiser_beta<F>(a: F) -> F
where
    F: Real + MulAdd<Output = F> + Pow<F, Output = F>,
    <F as Pow<F>>::Output: MulAdd<F, F>,
{
    if a > F::from(50).unwrap() {
        F::from(0.1102).unwrap() * (a - F::from(8.7).unwrap())
    } else if a > F::from(21).unwrap() {
        let a = a - F::from(21).unwrap();

        MulAdd::mul_add(
            a.pow(F::from(0.4).unwrap()),
            F::from(0.5842).unwrap(),
            F::from(0.07886).unwrap() * a,
        )
    } else if a > F::zero() {
        F::zero()
    } else {
        panic!("Expected a positive input.")
    }
}

/// Compute the attenuation of a Kaiser FIR filter.
///
/// Given the number of taps `N` and the transition width `width`, compute the
/// attenuation `a` in dB, given by Kaiser's formula:
/// ```custom
/// a = 2.285 * (N - 1) * pi * width + 7.95
/// ```
///
/// # Parameters
/// * `numtaps`: int  
///     The number of taps in the FIR filter.
/// * `width`: float  
///     The desired width of the transition region between passband and
///     stopband (or, in general, at any discontinuity) for the filter,
///     expressed as a fraction of the Nyquist frequency.
///
/// # Returns
/// `a`: The attenuation of the ripple, in dB.
///
/// # Examples
/// Suppose we want to design a FIR filter using the Kaiser window method
/// that will have 211 taps and a transition width of 9 Hz for a signal that
/// is sampled at 480 Hz. Expressed as a fraction of the Nyquist frequency,
/// the width is 9/(0.5*480) = 0.0375. The approximate attenuation (in dB)
/// is computed as follows:
/// ```
/// use sci_rs::signal::filter::design::kaiser_atten;
/// assert_eq!(64.48099630593983 , kaiser_atten(211, 0.0375));
/// ```
///
/// # See Also
/// [kaiserord], [kaiser_beta]
pub fn kaiser_atten<F>(numtaps: u32, width: F) -> F
where
    F: Real + MulAdd<Output = F> + RealField,
{
    MulAdd::mul_add(
        width,
        F::from(numtaps - 1).unwrap() * F::from(2.285).unwrap() * F::pi(),
        F::from(7.95).unwrap(),
    )
}

/// Determine the filter window parameters for the Kaiser window method.
///
/// The parameters returned by this function are generally used to create
/// a finite impulse response filter using the window method, with either
/// `firwin` or `firwin2`.
///
/// # Parameters
/// * `ripple`: float  
///     Upper bound for the deviation (in dB) of the magnitude of the
///     filter's frequency response from that of the desired filter (not
///     including frequencies in any transition intervals). That is, if w
///     is the frequency expressed as a fraction of the Nyquist frequency,
///     A(w) is the actual frequency response of the filter and D(w) is the
///     desired frequency response, the design requirement is that:
///     ```abs(A(w) - D(w))) < 10**(-ripple/20)```
///     for 0 <= w <= 1 and w not in a transition interval.
/// * `width`: float  
///     Width of transition region, normalized so that 1 corresponds to pi
///     radians / sample. That is, the frequency is expressed as a fraction
///     of the Nyquist frequency.
///
/// # Returns
/// * `numtaps`: int  
///     The length of the Kaiser window.
/// * `beta`: float  
///     The beta parameter for the Kaiser window.
///
/// # Notes
/// There are several ways to obtain the Kaiser window:
///
/// - ``signal.windows.kaiser(numtaps, beta, sym=True)``
/// - ``signal.get_window(beta, numtaps)``
/// - ``signal.get_window(('kaiser', beta), numtaps)``
///
/// The empirical equations discovered by Kaiser are used.
///
/// # References
/// Oppenheim, Schafer, "Discrete-Time Signal Processing", pp.475-476.
///
/// # Scipy Example
/// We will use the Kaiser window method to design a lowpass FIR filter
/// for a signal that is sampled at 1000 Hz.
///
/// We want at least 65 dB rejection in the stop band, and in the pass
/// band the gain should vary no more than 0.5%.
///
/// We want a cutoff frequency of 175 Hz, with a transition between the
/// pass band and the stop band of 24 Hz. That is, in the band [0, 163],
/// the gain varies no more than 0.5%, and in the band [187, 500], the
/// signal is attenuated by at least 65 dB.
///
/// ```custom,{class=language-python}
/// >>> import numpy as np
/// >>> from scipy.signal import kaiserord, firwin, freqz
/// >>> import matplotlib.pyplot as plt
/// >>> fs = 1000.0
/// >>> cutoff = 175
/// >>> width = 24
/// ```
///
/// The Kaiser method accepts just a single parameter to control the pass
/// band ripple and the stop band rejection, so we use the more restrictive
/// of the two. In this case, the pass band ripple is 0.005, or 46.02 dB,
/// so we will use 65 dB as the design parameter.
///
/// Use `kaiserord` to determine the length of the filter and the
/// parameter for the Kaiser window.
///
/// ```custom,{class=language-python}
/// >>> numtaps, beta = kaiserord(65, width/(0.5*fs))
/// >>> numtaps
/// 167
/// >>> beta
/// 6.20426
/// ```
///
/// Use `firwin` to create the FIR filter.
///
/// ```custom,{class=language-python}
/// >>> taps = firwin(numtaps, cutoff, window=('kaiser', beta),
/// ...               scale=False, fs=fs)
/// ```
///
/// Compute the frequency response of the filter.  ``w`` is the array of
/// frequencies, and ``h`` is the corresponding complex array of frequency
/// responses.
///
/// ```custom,{class=language-python}
/// >>> w, h = freqz(taps, worN=8000)
/// >>> w *= 0.5*fs/np.pi  # Convert w to Hz.
/// ```
///
/// Compute the deviation of the magnitude of the filter's response from
/// that of the ideal lowpass filter. Values in the transition region are
/// set to ``nan``, so they won't appear in the plot.
///
/// ```custom,{class=language-python}
/// >>> ideal = w < cutoff  # The "ideal" frequency response.
/// >>> deviation = np.abs(np.abs(h) - ideal)
/// >>> deviation[(w > cutoff - 0.5*width) & (w < cutoff + 0.5*width)] = np.nan
/// ```
///
/// Plot the deviation. A close look at the left end of the stop band shows
/// that the requirement for 65 dB attenuation is violated in the first lobe
/// by about 0.125 dB. This is not unusual for the Kaiser window method.
///
/// ```custom,{class=language-python}
/// >>> plt.plot(w, 20*np.log10(np.abs(deviation)))
/// >>> plt.xlim(0, 0.5*fs)
/// >>> plt.ylim(-90, -60)
/// >>> plt.grid(alpha=0.25)
/// >>> plt.axhline(-65, color='r', ls='--', alpha=0.3)
/// >>> plt.xlabel('Frequency (Hz)')
/// >>> plt.ylabel('Deviation from ideal (dB)')
/// >>> plt.title('Lowpass Filter Frequency Response')
/// >>> plt.show()
/// ```
///
/// See Also
/// --------
/// [kaiser_beta], [kaiser_atten]
///
pub fn kaiserord<F>(ripple: F, width: F) -> (F, F)
where
    F: Real + MulAdd<Output = F> + Pow<F, Output = F> + RealField,
{
    let a = Real::abs(ripple);
    if a < F::from(8).unwrap() {
        panic!("Requested maximum ripple attenuation is too small for the Kaiser formula.");
    }
    let beta = kaiser_beta(a);
    let numtaps =
        F::one() + (a - F::from(7.95).unwrap()) / (F::from(2.285).unwrap() * F::pi() * width);
    (Real::ceil(numtaps), beta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kaiserord_test_a() {
        let fs = 1000.;
        let cutoff = 175.;
        let width = 24.;
        let ripple = 65.;

        assert_eq!((167., 6.20426), kaiserord(ripple, width / (0.5 * fs)));
    }
}
