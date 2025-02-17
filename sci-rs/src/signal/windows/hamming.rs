use super::GeneralHamming;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Collection of arguments for window `Hamming` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct Hamming {
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl Hamming {
    /// Returns a Hamming struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `sym`: bool   
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, sym: bool) -> Self {
        Hamming { m, sym }
    }
}

#[cfg(feature = "alloc")]
impl<W> GetWindow<W> for Hamming
where
    W: Real + Float + RealField,
{
    /// Return a Hamming window.
    ///
    /// The Hamming window is a taper formed by using a raised cosine with non-zero endpoints,
    /// optimized to minimize the nearest side lobe.
    ///
    /// # Parameters
    /// * `M` : int  
    ///     Number of points in the output window. If zero, an empty array is returned. An
    ///     exception is thrown when it is negative.
    /// * `sym` : bool, optional  
    ///     When True (default), generates a symmetric window, for use in filter
    ///     design.  
    ///     When False, generates a periodic window, for use in spectral analysis.
    ///
    /// # Returns
    /// `w` : ndarray  
    ///     The window, with the maximum value normalized to 1 (though the value 1
    ///     does not appear if `M` is even and `sym` is True).
    ///
    /// # Notes
    /// The Hamming window is defined as
    ///
    /// $$w(n) = 0.54 - 0.46 \cos\left(\frac{2\pi{n}}{M-1}\right) \qquad 0 \leq n \leq M-1$$
    ///
    /// The Hamming was named for R. W. Hamming, an associate of J. W. Tukey and is described in
    /// Blackman and Tukey. It was recommended for smoothing the truncated autocovariance function
    /// in the time domain. Most references to the Hamming window come from the signal processing
    /// literature, where it is used as one of many windowing functions for smoothing values.  It
    /// is also known as an apodization (which means "removing the foot", i.e. smoothing
    /// discontinuities at the beginning and end of the sampled signal) or tapering function.
    ///
    /// # References
    /// [[1]] Blackman, R.B. and Tukey, J.W., (1958) The measurement of power spectra, Dover
    /// Publications, New York.  
    /// [[2]] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The University of Alberta
    /// Press, 1975, pp. 109-110.  
    /// [[3]] Wikipedia, "Window function", <https://en.wikipedia.org/wiki/Window_function>  
    /// [[4]] W.H. Press,  B.P. Flannery, S.A. Teukolsky, and W.T. Vetterling, "Numerical Recipes",
    /// Cambridge University Press, 1986, page 425.  
    /// [[5]] Scipy,
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.hamming.html>
    ///
    /// Examples
    /// --------
    /// Plot the window and its frequency response:
    ///
    /// ```custom,{class=language-python}
    /// >>> import numpy as np
    /// >>> from scipy import signal
    /// >>> from scipy.fft import fft, fftshift
    /// >>> import matplotlib.pyplot as plt
    /// ```
    ///
    /// ```custom,{class=language-python}
    /// >>> window = signal.windows.hamming(51)
    /// >>> plt.plot(window)
    /// >>> plt.title("Hamming window")
    /// >>> plt.ylabel("Amplitude")
    /// >>> plt.xlabel("Sample")
    /// ```
    ///
    /// ```custom,{class=language-python}
    /// >>> plt.figure()
    /// >>> A = fft(window, 2048) / (len(window)/2.0)
    /// >>> freq = np.linspace(-0.5, 0.5, len(A))
    /// >>> response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
    /// >>> plt.plot(freq, response)
    /// >>> plt.axis([-0.5, 0.5, -120, 0])
    /// >>> plt.title("Frequency response of the Hamming window")
    /// >>> plt.ylabel("Normalized magnitude [dB]")
    /// >>> plt.xlabel("Normalized frequency [cycles per sample]")
    /// ```
    ///
    /// The equivalent is:
    /// ```
    /// use sci_rs::signal::windows::{GetWindow, Hamming};
    /// let window: Vec<f64> = Hamming::new(51, true).get_window();
    /// ```
    ///
    /// [1]: #references
    /// [2]: #references
    /// [3]: #references
    /// [4]: #references
    /// [5]: #references
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        GeneralHamming::<W>::new(self.m, W::from(0.54).unwrap(), self.sym).get_window()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn hamming_case_a() {
        // from scipy.signal.windows import hamming
        // hamming(17)
        let h = Hamming::new(17, true);
        let expected = [
            0.08, 0.11501542, 0.21473088, 0.36396562, 0.54, 0.71603438, 0.86526912, 0.96498458, 1.,
            0.96498458, 0.86526912, 0.71603438, 0.54, 0.36396562, 0.21473088, 0.11501542, 0.08,
        ]
        .into();

        assert_vec_eq(expected, h.get_window());
    }

    #[test]
    fn hamming_case_b() {
        // from scipy.signal.windows import hamming
        // hamming(17, false)
        let h = Hamming::new(17, false);
        let expected = [
            0.08, 0.11106277, 0.2000559, 0.33496036, 0.49755655, 0.66588498, 0.81721193,
            0.93109988, 0.99216763, 0.99216763, 0.93109988, 0.81721193, 0.66588498, 0.49755655,
            0.33496036, 0.2000559, 0.11106277,
        ]
        .into();

        assert_vec_eq(expected, h.get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
