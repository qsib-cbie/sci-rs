use super::GeneralCosine;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Collection of arguments for window `Nuttall` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct Nuttall {
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl Nuttall {
    /// Returns a Nuttall struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `sym`: bool   
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, sym: bool) -> Self {
        Nuttall { m, sym }
    }
}

#[cfg(feature = "alloc")]
impl<W> GetWindow<W> for Nuttall
where
    W: Real + Float + RealField,
{
    /// Return a minimum 4-term Blackman-Harris window according to Nuttall.
    ///
    /// This variation is called "Nuttall4c" by Heinzel. [[2]]
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
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # References
    /// [[1]] A. Nuttall, "Some windows with very good sidelobe behavior," IEEE Transactions on
    /// Acoustics, Speech, and Signal Processing, vol. 29, no. 1, pp. 84-91, Feb 1981.
    /// :doi:`10.1109/TASSP.1981.1163506`.  
    /// [[2]] Heinzel G. et al., "Spectrum and spectral density estimation by the Discrete Fourier
    /// transform (DFT), including a comprehensive list of window functions and some new flat-top
    /// windows", February 15, 2002 <https://holometer.fnal.gov/GH_FFT.pdf>  
    /// [[3]] Scipy,
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.nuttall.html>
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
    /// >>> window = signal.windows.nuttall(51)
    /// >>> plt.plot(window)
    /// >>> plt.title("Nuttall window")
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
    /// use sci_rs::signal::windows::{GetWindow, Nuttall};
    /// let window: Vec<f64> = Nuttall::new(51, true).get_window();
    /// ```
    ///
    /// [1]: #references
    /// [2]: #references
    /// [3]: #references
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        GeneralCosine::<W>::new(
            self.m,
            [0.3635819, 0.4891775, 0.1365995, 0.0106411]
                .map(|f| W::from(f).unwrap())
                .into_iter()
                .collect(),
            self.sym,
        )
        .get_window()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn nuttall_case_a() {
        // from scipy.signal.windows import nuttall
        // nuttall(17)
        let h = Nuttall::new(17, true);
        let expected = [
            3.62800000e-04,
            4.15908007e-03,
            2.52055665e-02,
            8.96224370e-02,
            2.26982400e-01,
            4.44360497e-01,
            7.01958233e-01,
            9.16185585e-01,
            1.00000000e+00,
            9.16185585e-01,
            7.01958233e-01,
            4.44360497e-01,
            2.26982400e-01,
            8.96224370e-02,
            2.52055665e-02,
            4.15908007e-03,
            3.62800000e-04,
        ]
        .into();

        assert_vec_eq(expected, h.get_window());
    }

    #[test]
    fn nuttall_case_b() {
        // from scipy.signal.windows import nuttall
        // nuttall(17, false)
        let h = Nuttall::new(17, false);
        let expected = [
            3.62800000e-04,
            3.64256817e-03,
            2.10918726e-02,
            7.36770505e-02,
            1.87084736e-01,
            3.73448574e-01,
            6.11072447e-01,
            8.39394793e-01,
            9.80852709e-01,
            9.80852709e-01,
            8.39394793e-01,
            6.11072447e-01,
            3.73448574e-01,
            1.87084736e-01,
            7.36770505e-02,
            2.10918726e-02,
            3.64256817e-03,
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
