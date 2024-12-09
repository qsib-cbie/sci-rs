use super::GeneralCosine;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Collection of arguments for window `GeneralHamming` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct GeneralHamming<F>
where
    F: Real,
{
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// The window coefficient, α.
    pub alpha: F,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl<F> GeneralHamming<F>
where
    F: Real,
{
    /// Returns a GeneralHamming struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `alpha` : float  
    ///     The window coefficient, α.
    /// * `sym`: bool   
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, alpha: F, sym: bool) -> Self {
        GeneralHamming { m, alpha, sym }
    }
}

#[cfg(feature = "alloc")]
impl<F, W> GetWindow<W> for GeneralHamming<F>
where
    F: Real,
    W: Real + Float + RealField,
{
    /// Return a generalized Hamming window.
    ///
    /// The generalized Hamming window is constructed by multiplying a rectangular
    /// window by one period of a cosine function [[1][1]].
    ///
    /// # Parameters
    /// * `M` : usize  
    ///     Number of points in the output window. If zero, an empty array is returned. An
    ///     exception is thrown when it is negative.
    /// * `alpha` : float  
    ///     The window coefficient, α.
    /// * `sym` : bool, optional  
    ///     When True (default), generates a symmetric window, for use in filter design.  
    ///     When False, generates a periodic window, for use in spectral analysis.
    ///
    /// # Returns
    /// `w` : ndarray  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # See Also
    /// hamming, hann
    ///
    /// # Notes
    /// The generalized Hamming window is defined as  
    /// $$w(n) = \alpha - \left(1 - \alpha\right)
    ///           \cos\left(\frac{2\pi{n}}{M-1}\right) \qquad 0 \leq n \leq M-1$$  
    /// Both the common Hamming window and Hann window are special cases of the generalized Hamming
    /// window with :math:`\alpha` = 0.54 and :math:`\alpha` = 0.5, respectively [[2][2]].
    ///
    /// # References
    ///
    /// [[1]] DSPRelated, "Generalized Hamming Window Family",
    /// <https://www.dsprelated.com/freebooks/sasp/Generalized_Hamming_Window_Family.html>  
    /// [[2]] Wikipedia, "Window function", <https://en.wikipedia.org/wiki/Window_function>  
    /// [[3]] Riccardo Piantanida ESA, "Sentinel-1 Level 1 Detailed Algorithm Definition",
    /// <https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Level-1-Detailed-Algorithm-Definition>  
    /// [[4]] Matthieu Bourbigot ESA, "Sentinel-1 Product Definition",
    /// <https://sentinel.esa.int/documents/247904/1877131/Sentinel-1-Product-Definition>  
    /// [[5]] Scipy,
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.general_hamming.html>
    ///
    /// # Examples
    /// The Sentinel-1A/B Instrument Processing Facility uses generalized Hamming
    /// windows in the processing of spaceborne Synthetic Aperture Radar (SAR)
    /// data [[3][3]]. The facility uses various values for the :math:`\alpha`
    /// parameter based on operating mode of the SAR instrument. Some common
    /// :math:`\alpha` values include 0.75, 0.7 and 0.52 [[4][4]]. As an example, we
    /// plot these different windows.
    ///
    /// ```custom,{class=language-python}
    /// >>> import numpy as np
    /// >>> from scipy.signal.windows import general_hamming
    /// >>> from scipy.fft import fft, fftshift
    /// >>> import matplotlib.pyplot as plt
    /// ```
    ///
    /// ```custom,{class=language-python}
    /// >>> fig1, spatial_plot = plt.subplots()
    /// >>> spatial_plot.set_title("Generalized Hamming Windows")
    /// >>> spatial_plot.set_ylabel("Amplitude")
    /// >>> spatial_plot.set_xlabel("Sample")
    /// ```
    ///
    /// ```custom,{class=language-python}
    /// >>> fig2, freq_plot = plt.subplots()
    /// >>> freq_plot.set_title("Frequency Responses")
    /// >>> freq_plot.set_ylabel("Normalized magnitude [dB]")
    /// >>> freq_plot.set_xlabel("Normalized frequency [cycles per sample]")
    /// ```
    ///
    /// ```custom,{class=language-python}
    /// >>> for alpha in [0.75, 0.7, 0.52]:
    /// ...     window = general_hamming(41, alpha)
    /// ...     spatial_plot.plot(window, label="{:.2f}".format(alpha))
    /// ...     A = fft(window, 2048) / (len(window)/2.0)
    /// ...     freq = np.linspace(-0.5, 0.5, len(A))
    /// ...     response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
    /// ...     freq_plot.plot(freq, response, label="{:.2f}".format(alpha))
    /// >>> freq_plot.legend(loc="upper right")
    /// >>> spatial_plot.legend(loc="upper right")
    /// ```
    ///
    /// The equivalent is:
    /// ```
    /// use sci_rs::signal::windows::{GetWindow, GeneralHamming};
    /// let window = GeneralHamming::new(41, 0.75, true);
    /// ```
    ///
    /// [1]: #references
    /// [2]: #references
    /// [3]: #references
    /// [4]: #references
    /// [5]: #references
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        GeneralCosine::new(self.m, [self.alpha, F::one() - self.alpha].into(), self.sym)
            .get_window()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn general_hamming_case_a() {
        // from scipy.signal.windows import general_hamming
        // general_hamming(15, 0.65)
        let gh = GeneralHamming::new(15, 0.65, true);

        let expected = [
            0.3, 0.3346609, 0.43177857, 0.57211767, 0.72788233, 0.86822143, 0.9653391, 1.,
            0.9653391, 0.86822143, 0.72788233, 0.57211767, 0.43177857, 0.3346609, 0.3,
        ]
        .into();
        assert_vec_eq(expected, gh.get_window());
    }

    #[test]
    fn general_hamming_case_b() {
        // from scipy.signal.windows import general_hamming
        // general_hamming(15, 0.65)
        let gh = GeneralHamming::new(15, 0.65, false);

        let expected = [
            0.3, 0.33025909, 0.41580429, 0.54184405, 0.68658496, 0.825, 0.93315595, 0.99235166,
            0.99235166, 0.93315595, 0.825, 0.68658496, 0.54184405, 0.41580429, 0.33025909,
        ]
        .into();
        assert_vec_eq(expected, gh.get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
