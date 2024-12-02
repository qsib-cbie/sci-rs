use super::{extend, len_guard, truncate};
use crate::special::Bessel;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;

#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `Kaiser` for use in [GetWindow].
#[cfg(feature = "alloc")]
#[derive(Debug, Clone, PartialEq)]
pub struct Kaiser<F>
where
    F: Real,
{
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Shape parameter.
    pub beta: F,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl<F> Kaiser<F>
where
    F: Real,
{
    /// Returns a Kaiser struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `beta` : float  
    ///     Shape parameter, determines trade-off between main-lobe width and side lobe level. As
    ///     beta gets large, the window narrows.
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, beta: F, sym: bool) -> Self {
        Kaiser { m, beta, sym }
    }
}

impl<F, W> GetWindow<W> for Kaiser<F>
where
    F: Real,
    W: Real + Bessel,
{
    /// Return a [Kaiser] window.
    ///
    /// The Kaiser window is a taper formed by using a Bessel function.
    ///
    /// Parameters
    /// ----------
    /// * `M` : int  
    ///     Number of points in the output window. If zero, an empty array is returned. An
    ///     exception is thrown when it is negative.
    /// * `beta` : float  
    ///     Shape parameter, determines trade-off between main-lobe width and side lobe level. As
    ///     beta gets large, the window narrows.
    /// * `sym` : bool, optional  
    ///     When True (default), generates a symmetric window, for use in filter design.  
    ///     When False, generates a periodic window, for use in spectral analysis.
    ///
    /// Returns
    /// -------
    /// `w` : ndarray  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// Notes
    /// -----
    /// The Kaiser window is defined as
    ///
    /// $$w(n) = I_0\left( \beta \sqrt{1-\frac{4n^2}{(M-1)^2}} \right)/I_0(\beta)$$
    ///
    /// with
    ///
    /// $$\quad -\frac{M-1}{2} \leq n \leq \frac{M-1}{2}$$,
    ///
    /// where $I_0$ is the modified zeroth-order Bessel function.
    ///
    /// The Kaiser was named for Jim Kaiser, who discovered a simple approximation to the DPSS
    /// window based on Bessel functions.  
    /// The Kaiser window is a very good approximation to the Digital Prolate Spheroidal Sequence,
    /// or Slepian window, which is the transform which maximizes the energy in the main lobe of
    /// the window relative to total energy.
    ///
    /// The Kaiser can approximate other windows by varying the beta parameter. [(Some literature
    /// uses alpha = beta/pi.)]([4])
    ///
    /// ====  =======================  
    /// beta  Window shape  
    /// ====  =======================  
    /// 0     Rectangular  
    /// 5     Similar to a Hamming  
    /// 6     Similar to a Hann  
    /// 8.6   Similar to a Blackman  
    /// ====  =======================  
    ///
    /// A beta value of 14 is probably a good starting point. Note that as beta gets large, the
    /// window narrows, and so the number of samples needs to be large enough to sample the
    /// increasingly narrow spike, otherwise NaNs will be returned.
    ///
    /// Most references to the Kaiser window come from the signal processing literature, where it
    /// is used as one of many windowing functions for smoothing values.  It is also known as an
    /// apodization (which means "removing the foot", i.e. smoothing discontinuities at the
    /// beginning and end of the sampled signal) or tapering function.
    ///
    /// # References
    /// [[1]] J. F. Kaiser, "Digital Filters" - Ch 7 in "Systems analysis by digital computer",
    /// Editors: F.F. Kuo and J.F. Kaiser, p 218-285. John Wiley and Sons, New York, (1966).  
    /// [[2]] E.R. Kanasewich, "Time Sequence Analysis in Geophysics", The University of Alberta
    /// Press, 1975, pp. 177-178.  
    /// [[3]] Wikipedia, "Window function", <https://en.wikipedia.org/wiki/Window_function>  
    /// [[4]] F. J. Harris, "On the use of windows for harmonic analysis with the discrete Fourier
    /// transform," Proceedings of the IEEE, vol. 66, no. 1, pp. 51-83, Jan. 1978.  
    /// :doi:`10.1109/PROC.1978.10837`.
    /// [[5]]
    /// [Scipy](<https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.kaiser.html>)
    ///
    /// # Examples
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
    /// >>> window = signal.windows.kaiser(51, beta=14)
    /// >>> plt.plot(window)
    /// >>> plt.title(r"Kaiser window ($\beta$=14)")
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
    /// >>> plt.title(r"Frequency response of the Kaiser window ($\beta$=14)")
    /// >>> plt.ylabel("Normalized magnitude [dB]")
    /// >>> plt.xlabel("Normalized frequency [cycles per sample]")
    /// ```
    ///
    /// The equivalent is:
    /// ```
    /// use sci_rs::signal::windows::{GetWindow, Kaiser};
    /// let window: Vec<f64> = Kaiser::new(51, 14., true).get_window();
    /// ```
    ///
    /// Due to current implementation limitations, note that the following might not possible:
    /// ```
    /// use sci_rs::signal::windows::{GetWindow, Kaiser};
    /// let window: Vec<f32> = Kaiser::new(51, 14., true).get_window();
    /// println!("window = {:?}", window);
    /// ```
    ///
    /// [1]: #references
    /// [2]: #references
    /// [3]: #references
    /// [4]: #references
    /// [5]: #references
    fn get_window(&self) -> Vec<W> {
        if len_guard(self.m) {
            return Vec::<W>::new();
        }
        let (m, needs_trunc) = extend(self.m, self.sym);
        let n = (0..m);
        let alpha = W::from(self.m - 1).unwrap() / W::from(2).unwrap();
        let beta = W::from(self.beta).unwrap();
        // let w: Vec<W> = n
        //     .map(|ni| W::from(ni).unwrap() - alpha)
        //     .map(|ni| (ni / alpha).powf(W::from(2).unwrap()))
        //     .map(|ni| (W::one() - ni).sqrt())
        //     .map(|ni| i0(beta * ni) / i0(beta))
        //     .collect();
        let w: Vec<W> = n
            .map(|ni| {
                (beta
                    * (W::one()
                        - ((W::from(ni).unwrap() - alpha) / alpha).powf(W::from(2).unwrap()))
                    .sqrt())
                .i0()
                    / (beta.i0())
            })
            .collect();
        truncate(w, needs_trunc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn kaiser_17_8_true() {
        // from scipy.signal.windows import kaiser
        // kaiser(17, beta = 0.8)
        let expected = vec![
            0.85725436, 0.88970403, 0.9183205, 0.94289618, 0.96325245, 0.97924114, 0.99074569,
            0.99768219, 1., 0.99768219, 0.99074569, 0.97924114, 0.96325245, 0.94289618, 0.9183205,
            0.88970403, 0.85725436,
        ];
        let k = Kaiser::new(17, 0.8, true);

        assert_vec_eq(expected, k.get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
