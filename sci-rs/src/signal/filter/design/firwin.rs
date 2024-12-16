use super::filter_type::FilterBandType;
use super::iirfilter_dyn;
use super::{kaiser_atten, kaiser_beta};
use crate::signal::windows::get_window;
use crate::{error, error::Error, special};
use core::cmp::Ord;
use nalgebra::RealField;
use num_traits::{real::Real, Float, MulAdd, Pow};

#[cfg(feature = "alloc")]
use crate::signal::{
    windows,
    windows::{GetWindow, GetWindowBuilder},
    windows::{Hamming, Kaiser},
};
#[cfg(feature = "alloc")]
use alloc::{string::String, vec, vec::Vec};

/// Validation of [firwin_dyn] input.
fn firwin_dyn_validate<F: Real + PartialOrd, W: Real>(
    numtaps: &usize,
    cutoff: &[F],
    pass_zero: &FilterBandType,
    width: &Option<F>,
    window: &Option<&impl GetWindow<W>>,
) -> Result<(), Error> {
    if cutoff.is_empty() {
        return Err(Error::InvalidArg {
            arg: "cutoff".into(),
            reason: "At least one cutoff frequency must be given.".into(),
        });
    }

    // Whilst it may be faster to write
    // if *cutoff.iter().min().unwrap() <= F::zero() || *cutoff.iter().max().unwrap() >= F::one() {
    //
    // vec<f64>.min() requires the Ord trait, which is very difficult to impose onto f64
    let minimal = cutoff.iter().fold(F::max_value(), |acc, &m| acc.min(m));
    let maximal = cutoff.iter().fold(F::min_value(), |acc, &m| acc.max(m));
    if minimal <= F::zero() || maximal >= F::one() {
        return Err(Error::InvalidArg {
            arg: "cutoff".into(),
            reason:
                "Invalid cutoff frequency: frequencies must be greater than 0 and less than fs/2."
                    .into(),
        });
    }
    if cutoff.windows(2).any(|x| x[0] >= x[1]) {
        return Err(Error::InvalidArg {
            arg: "cutoff".into(),
            reason: "Invalid cutoff frequencies: the frequencies must be strictly increasing."
                .into(),
        });
    }

    match pass_zero {
        FilterBandType::Lowpass => {
            if cutoff.len() != 1 {
                return Err(Error::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have one element if pass_zero is Lowpass.".into(),
                });
            }
        }
        FilterBandType::Bandstop => {
            if cutoff.len() < 2 {
                return Err(Error::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have at least two elements if pass_zero is Bandstop."
                        .into(),
                });
            }
        }
        FilterBandType::Highpass => {
            if cutoff.len() != 1 {
                return Err(Error::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have one element if pass_zero is Highpass.".into(),
                });
            }
        }
        FilterBandType::Bandpass => {
            if cutoff.len() < 2 {
                return Err(Error::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have at least two elements if pass_zero is Bandpass."
                        .into(),
                });
            }
        }
    }

    // While this was silently ignored in Scipy, we make this explicit here.
    // Cannot use != on &impl GetWindow.
    if window.is_some() && width.is_some() {
        return Err(Error::InvalidArg {
            arg: "cutoff".into(),
            reason: "Setting both window and with to something is silently ignored only in Scipy."
                .into(),
        });
    }

    // This is here only because impl GetWindow on GetWindowBuilder is non-trivial to fetch numtaps
    // outside the struct/enum-variants.
    // if let Some(w) = *window {
    //     if w.m != numtaps {
    //         return Err(Error::ConfictArg {
    //             arg: "window".into(),
    //             reason: "Window has m value differing from numtaps",
    //         });
    //     }
    // }
    Ok(())
}

/// FIR filter design using the window method.
///
/// This function computes the coefficients of a finite impulse response filter. The filter will
/// have linear phase; it will be Type I if `numtaps` is odd and Type II if `numtaps` is even.
///
/// Type II filters always have zero response at the Nyquist frequency, so a [error::Error] is
/// returned if firwin is called with `numtaps` even and having a passband whose right end is at
/// the Nyquist frequency.
///
/// # Parameters
/// * `numtaps`: usize  
///     Length of the filter (number of coefficients, i.e. the filter order + 1).  `numtaps` must
///     be odd if a passband includes the Nyquist frequency.
/// * `cutoff`: 1-D array_like  
///     Cutoff frequency of filter (expressed in the same units as `fs`) OR an array of cutoff
///     frequencies (that is, band edges). In the latter case, the frequencies in `cutoff` should
///     be positive and monotonically increasing between 0 and `fs/2`. The values 0 and `fs/2` must
///     not be included in `cutoff`.
/// * `width`: float or None, optional  
///     If `width` is not None, then assume it is the approximate width of the transition region
///     (expressed in the same units as `fs`) for use in Kaiser FIR filter design. In this case,
///     the `window` argument is ignored.
/// * `window` : string or tuple of string and parameter values, optional  
///     Desired window to use. See [GetWindow] for a list of windows and required parameters.  
///     Defaults to Hamming.
///     Please set `sym=True` if you wish to have similar behaviour to scipy.
/// * `pass_zero` : [FilterBandType], optional  
///     If True, the gain at the frequency 0 (i.e., the "DC gain") is 1.  
///     If False, the DC gain is
///     0. Can also be a string argument for the desired filter type (equivalent to ``btype`` in
///        [IIR design functions](iirfilter_dyn)).
///
/// * `scale` : bool, optional  
///     Set to True to scale the coefficients so that the frequency response is exactly unity at a
///     certain frequency. That frequency is either:
///
///     * 0 (DC) if the first passband starts at 0 (i.e. pass_zero
///       is True)
///     * `fs/2` (the Nyquist frequency) if the first passband ends at
///       `fs/2` (i.e the filter is a single band highpass filter);
///       center of first passband otherwise
///
/// * fs : float, optional  
///     The sampling frequency of the signal. Each frequency in `cutoff`
///     must be between 0 and ``fs/2``.  Default is 2.
///
/// Returns
/// -------
/// `h` : Result<`W`, _>
///     Coefficients of length `numtaps` FIR filter.
///
/// Raises
/// ------
/// Error
///     If any value in `cutoff` is less than or equal to 0 or greater than or equal to ``fs/2``,
///     if the values in `cutoff` are not strictly monotonically increasing, or if `numtaps` is
///     even but a passband includes the Nyquist frequency.
///
/// See Also
/// --------
/// firwin2
/// firls
/// minimum_phase
/// remez
///
/// Examples
/// --------
/// * Low-pass from 0 to f:
///
/// ```custom,{class=language-python}
/// >>> from scipy import signal
/// >>> numtaps = 3
/// >>> f = 0.1
/// >>> signal.firwin(numtaps, f)
/// array([ 0.06799017,  0.86401967,  0.06799017])
/// ```
/// Sci-rs:
/// ```
/// use approx:: assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Lowpass,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![0.06799017,  0.86401967,  0.06799017];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// Use a specific window function:
///
/// ```custom,{class=language-python}
/// >>> signal.firwin(numtaps, f, window='nuttall')
/// array([  3.56607041e-04,   9.99286786e-01,   3.56607041e-04])
/// ```
/// Sci-rs:
//  TODO: This still needs work so that its not Some(&Nuttall::new(...)) which might be mistakenly
//  different from what's above it.
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Nuttall;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64],
///     None,
///     Some(&Nuttall::new(3, true)),
///     &FilterBandType::Lowpass,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![3.56607041e-04, 9.99286786e-01, 3.56607041e-04];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// High-pass ('stop' from 0 to f):
///
/// ```custom,{class=language-python}
/// >>> signal.firwin(numtaps, f, pass_zero=False)
/// array([-0.00859313,  0.98281375, -0.00859313])
/// ```
/// Sci-rs:
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Highpass,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![-0.00859313, 0.98281375, -0.00859313];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// Band-pass:
///
/// ```custom,{class=language-python}
/// >>> f1, f2 = 0.1, 0.2
/// >>> signal.firwin(numtaps, [f1, f2], pass_zero=False)
/// array([ 0.06301614,  0.88770441,  0.06301614])
/// ```
/// Sci-rs:
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64, 0.2f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Bandpass,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![ 0.06301614,  0.88770441,  0.06301614];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// Band-stop:
///
/// ```custom,{class=language-python}
/// >>> signal.firwin(numtaps, [f1, f2])
/// array([-0.00801395,  1.0160279 , -0.00801395])
/// ```
/// Sci-rs:
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64, 0.2f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Bandstop,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![-0.00801395, 1.0160279 , -0.00801395];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// Multi-band (passbands are `[0, f1]`, `[f2, f3]` and `[f4, 1]`):
///
/// ```custom,{class=language-python}
/// >>> f3, f4 = 0.3, 0.4
/// >>> signal.firwin(numtaps, [f1, f2, f3, f4])
/// array([-0.01376344,  1.02752689, -0.01376344])
/// ```
/// Sci-rs:
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64, 0.2f64, 0.3f64, 0.4f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Bandstop,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![-0.01376344, 1.02752689, -0.01376344];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
///
/// Multi-band (passbands are `[f1, f2]` and `[f3,f4]`):
///
/// ```custom,{class=language-python}
/// >>> signal.firwin(numtaps, [f1, f2, f3, f4], pass_zero=False)
/// array([ 0.04890915,  0.91284326,  0.04890915])
/// ```
/// Sci-rs:
/// ```
/// use approx::assert_abs_diff_eq;
/// use sci_rs::signal::filter::design::{firwin_dyn, FilterBandType};
/// use sci_rs::signal::windows::Hamming;
///
/// let window: Vec<f64> = firwin_dyn(
///     3,
///     &[0.1f64, 0.2f64, 0.3f64, 0.4f64],
///     None,
///     None::<&Hamming>,
///     &FilterBandType::Bandpass,
///     None,
///     None,
/// )
/// .unwrap();
/// let expected = vec![0.04890915, 0.91284326, 0.04890915];
///
/// fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
///     for (a, b) in a.into_iter().zip(b) {
///         assert_abs_diff_eq!(a, b, epsilon = 1e-6);
///     }
/// }
///
/// assert_vec_eq(window, expected);
/// ```
// In accordance with https://github.com/scipy/scipy/pull/16315, nyq as an argument should not be
// provided.
#[cfg(feature = "alloc")]
pub fn firwin_dyn<F, W>(
    numtaps: usize,
    cutoff: &[F],
    width: Option<F>,
    window: Option<&impl GetWindow<W>>,
    pass_zero: &FilterBandType, // Union with bools to follow scipy?
    scale: Option<bool>,
    fs: Option<F>,
) -> Result<Vec<W>, Error>
where
    W: Real + Float + RealField + special::Bessel,
    F: Real + PartialOrd + Float + RealField + MulAdd<Output = F> + Pow<F, Output = F>,
{
    let fs = match fs {
        None => F::from(2).unwrap(),
        Some(x) => x,
    };
    let nyq = fs / F::from(2).unwrap();
    let cutoff: Vec<_> = cutoff.iter().map(|&c| c / nyq).collect();

    firwin_dyn_validate(&numtaps, &cutoff, pass_zero, &width, &window)?;

    // Get and apply the window function.
    let win: Vec<W> = if let Some(x) = window {
        // ?warn if window.sym != true
        x.get_window()
    } else if let Some(width) = width {
        let atten = kaiser_atten(numtaps.try_into().unwrap(), width / nyq);
        let beta = kaiser_beta(atten);
        let k = get_window(GetWindowBuilder::Kaiser { beta }, numtaps, Some(false));
        k.get_window()
    } else {
        let h = get_window(GetWindowBuilder::<F>::Hamming, numtaps, Some(false));
        h.get_window()
    };

    let pass_zero: bool = match pass_zero {
        FilterBandType::Lowpass => true,
        FilterBandType::Bandstop => true,
        FilterBandType::Highpass => false,
        FilterBandType::Bandpass => false,
    };
    let pass_nyquist = (cutoff.len() % 2 == 1) ^ pass_zero;
    if pass_nyquist && numtaps % 2 == 0 {
        return Err(Error::InvalidArg {
                    arg: "numtaps".into(),
                    reason: "A filter with an even number of coefficients must have zero response at the Nyquist frequency."
                        .into(),
                });
    }

    let cutoff: Vec<F> = {
        let mut tmp = Vec::<F>::new();
        if pass_zero {
            tmp.push(F::zero());
        }
        tmp.extend_from_slice(&cutoff);
        if pass_nyquist {
            tmp.push(F::one());
        }
        tmp
    };
    if cutoff.len() % 2 != 0 {
        unreachable!();
        // return Err(Error::InvalidArg {
        //     arg: "cutoff".into(),
        //     reason: "Parity of cutoff given for type of Filter is incorrect.".into(),
        // });
    }
    let bands: Vec<_> = cutoff.chunks_exact(2).collect();
    let scale_frequency = {
        // for use in scale branch later
        let left = bands[0][0];
        let right = bands[0][1];
        if left == F::zero() {
            F::zero()
        } else if right == F::one() {
            F::one()
        } else {
            F::from(0.5).unwrap() * (left + right)
        }
    };

    // Build up the coefficients.
    // alpha: scalar = 0.5 * (numtaps - 1)
    let alpha = F::from(0.5).unwrap() * F::from(numtaps - 1).unwrap();
    // m: 1Darray[numtaps] = 0..numtaps - alpha // lifetimes
    // h: 1Darray[numtaps] = sum([right * sinc(right *m) - left * sinc(left * m)) for (left, right)
    //     in bands])
    let h: Vec<W> = bands
        .into_iter()
        .map(|b| {
            let left = b[0];
            let right = b[1];
            let m = (0..numtaps).map(|mi| F::from(mi).unwrap() - alpha);
            let h: Vec<F> = m
                .map(|mi| {
                    if !mi.is_zero() {
                        right * (right * mi * F::pi()).sinc() - left * (left * mi * F::pi()).sinc()
                    } else {
                        right - left
                    }
                })
                .collect();
            h
        })
        .fold(vec![W::zero(); numtaps], |mut acc, x| {
            acc.iter_mut()
                .zip(x)
                .map(|(a, xi)| *a += W::from(xi).unwrap())
                .collect::<()>();
            acc
        });
    let mut h: Vec<W> = h.into_iter().zip(win).map(|(hi, wi)| hi * wi).collect();

    let scale = scale.unwrap_or(true);
    if scale {
        // m: 1Darray[numtaps] = 0..numtaps - alpha // lifetimes
        let m = (0..numtaps).map(|mi| F::from(mi).unwrap() - alpha);
        // c: 1Darray[numtaps] = np.cos(np.pi * m * scale_frequency)
        let c = m
            .map(|mi| F::pi() * mi * scale_frequency)
            .map(|x| W::from(Real::cos(x)).unwrap());
        // s: scalar = np.sum(h*c)
        let s: W = h
            .iter()
            .zip(c)
            .map(|(&hi, ci)| hi * ci)
            .fold(W::zero(), |acc, x| acc + x);
        // h: 1Darray[numtaps] = h/s
        h.iter_mut().map(|hi| *hi /= s).collect::<()>();
    }

    Ok(h)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::signal::windows;
    use approx::assert_abs_diff_eq;

    #[test]
    fn invalid_args() {
        type E = crate::error::Error;
        {
            let cutoff = Vec::<f64>::new();
            let empty_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Lowpass,
                None,
                None,
            );
            assert_eq!(
                empty_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "At least one cutoff frequency must be given.".into()
                }
            );
        }
        {
            let cutoff: Vec<f64> = vec![-0.1, 3.];
            let invalid_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Bandpass,
                None,
                None,
            );
            assert_eq!(
                invalid_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "Invalid cutoff frequency: frequencies must be greater than 0 and less than fs/2.".into()
                }
            );
        }
        {
            let cutoff: Vec<f64> = vec![0.2, 0.2];
            let decreasing_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Bandpass,
                None,
                None,
            );
            assert_eq!(
                decreasing_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason:
                        "Invalid cutoff frequencies: the frequencies must be strictly increasing."
                            .into()
                }
            );
        }
        {
            // Lowpass
            let cutoff: Vec<f64> = vec![0.2, 0.7];
            let lowpass_invalid_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Lowpass,
                None,
                None,
            );
            assert_eq!(
                lowpass_invalid_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have one element if pass_zero is Lowpass.".into()
                }
            );
        }
        {
            // Bandstop
            let cutoff: Vec<f64> = vec![0.7];
            let bandstop_invalid_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Bandstop,
                None,
                None,
            );
            assert_eq!(
                bandstop_invalid_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have at least two elements if pass_zero is Bandstop."
                        .into()
                }
            );
        }
        {
            // Highpass
            let cutoff: Vec<f64> = vec![0.2, 0.7];
            let highpass_invalid_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Highpass,
                None,
                None,
            );
            assert_eq!(
                highpass_invalid_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have one element if pass_zero is Highpass.".into()
                }
            );
        }
        {
            // Bandpass
            let cutoff: Vec<f64> = vec![0.2];
            let bandpass_invalid_cutoff: Result<Vec<f64>, E> = firwin_dyn(
                3,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Bandpass,
                None,
                None,
            );
            assert_eq!(
                bandpass_invalid_cutoff.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "cutoff must have at least two elements if pass_zero is Bandpass."
                        .into()
                }
            );
        }
        {
            let cutoff: Vec<f64> = vec![0.2, 0.7];
            let invalid_numtaps: Result<Vec<f64>, E> = firwin_dyn(
                4,
                &cutoff,
                None,
                None::<&Hamming>,
                &FilterBandType::Bandstop,
                None,
                None,
            );
            assert_eq!(
                invalid_numtaps.unwrap_err(),
                E::InvalidArg {
                    arg: "numtaps".into(),
                    reason: "A filter with an even number of coefficients must have zero response at the Nyquist frequency.".into()
                }
            );
        }
    }

    #[test]
    fn conflicting_args() {
        type E = crate::error::Error;
        {
            let cutoff: Vec<f64> = vec![0.2, 0.5, 0.7];
            let window = Hamming::new(12, true);
            let conflicting: Result<Vec<f64>, E> = firwin_dyn(
                5,
                &cutoff,
                Some(0.3),
                Some(&window),
                &FilterBandType::Bandpass,
                None,
                None,
            );
            assert_eq!(
                conflicting.unwrap_err(),
                E::InvalidArg {
                    arg: "cutoff".into(),
                    reason: "Setting both window and with to something is silently ignored only in Scipy."
                        .into()
                }
            );
        }
    }

    #[test]
    fn bandstop() {
        // from scipy.signal import firwin
        // firwin(numtaps=5, cutoff=[0.2, 0.7], width=None, # window=None,
        //        pass_zero='bandstop', scale=True, fs=None)
        let cutoff: Vec<f64> = vec![0.2, 0.7];
        let window: Result<Vec<f64>, _> = firwin_dyn(
            5,
            &cutoff,
            None,
            None::<&Hamming>,
            &FilterBandType::Bandstop,
            None,
            None,
        );
        let expected = vec![0.05126868, -0.08050021, 1.05846306, -0.08050021, 0.05126868];
        assert_vec_eq(expected, window.unwrap());
    }

    #[test]
    fn bandpass() {
        // from scipy.signal import firwin
        // firwin(numtaps=5, cutoff=[0.2, 0.7], width=None, # window=None,
        //        pass_zero='bandpass', scale=True, fs=None)
        let cutoff: Vec<f64> = vec![0.2, 0.7];
        let window: Result<Vec<f64>, _> = firwin_dyn(
            5,
            &cutoff,
            None,
            None::<&Hamming>,
            &FilterBandType::Bandpass,
            None,
            None,
        );
        let expected = vec![-0.04340507, 0.06815306, 0.89611567, 0.06815306, -0.04340507];
        assert_vec_eq(expected, window.unwrap());
    }

    #[test]
    fn lowpass() {
        // from scipy.signal import firwin
        // firwin(numtaps=5, cutoff=[0.2], width=None, # window=None,
        //        pass_zero='lowpass', scale=True, fs=None)
        let cutoff: Vec<f64> = vec![0.2];
        let window: Result<Vec<f64>, _> = firwin_dyn(
            5,
            &cutoff,
            None,
            None::<&Hamming>,
            &FilterBandType::Lowpass,
            None,
            None,
        );
        let expected = vec![0.02840647, 0.23700821, 0.46917063, 0.23700821, 0.02840647];
        assert_vec_eq(expected, window.unwrap());
    }

    #[test]
    fn highpass() {
        // from scipy.signal import firwin
        // firwin(numtaps=5, cutoff=[0.2], width=None, # window=None,
        //        pass_zero='highpass', scale=True, fs=None)
        let cutoff: Vec<f64> = vec![0.2];
        let window: Result<Vec<f64>, _> = firwin_dyn(
            5,
            &cutoff,
            None,
            None::<&Hamming>,
            &FilterBandType::Highpass,
            None,
            None,
        );
        let expected = vec![-0.01238356, -0.1033217, 0.81812371, -0.1033217, -0.01238356];
        assert_vec_eq(expected, window.unwrap());
    }

    #[test]
    fn different_fs() {
        // from scipy.signal import firwin
        // firwin(numtaps=5, cutoff=[2], width=None, # window=None,
        //        pass_zero='lowpass', scale=True, fs=10)
        let cutoff: Vec<f64> = vec![2.];
        let window: Result<Vec<f64>, _> = firwin_dyn(
            5,
            &cutoff,
            None,
            None::<&Hamming>,
            &FilterBandType::Lowpass,
            None,
            Some(10.),
        );
        let expected = vec![0.01008727, 0.22034079, 0.53914388, 0.22034079, 0.01008727];
        assert_vec_eq(expected, window.unwrap());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
