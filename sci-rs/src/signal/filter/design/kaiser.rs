use core::f64::consts::PI;
use core::ops::Mul;
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
    } else if a > F::from(0).unwrap() {
        F::from(0).unwrap()
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
    F: Real + MulAdd<Output = F>,
{
    MulAdd::mul_add(
        width,
        F::from(numtaps - 1).unwrap() * F::from(2.285 * PI).unwrap(),
        F::from(7.95).unwrap(),
    )
}
