use super::GeneralCosine;
use nalgebra::RealField;
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `Blackman` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct Blackman {
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl Blackman {
    /// Returns a Blackman struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, sym: bool) -> Self {
        Blackman { m, sym }
    }
}

#[cfg(feature = "alloc")]
impl<W> GetWindow<W> for Blackman
where
    W: Real + Float + RealField,
{
    /// Return a window of type: Blackman.
    ///
    /// The Blackman window is a taper formed by using the first three terms of a summation of
    /// cosines. It was designed to have close to the minimal leakage possible. It is close to
    /// optimal, only slightly worse than a Kaiser window.
    ///
    /// # Parameters
    /// `self`: [Blackman]
    ///
    /// # Returns
    /// `w`: `vec<F>`  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # Example
    /// ```
    /// use sci_rs::signal::windows::{Blackman, GetWindow};
    ///
    /// let nx = 8;
    /// let b = Blackman::new(nx, true);
    /// ```
    ///
    /// # References
    /// <https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.blackman.html>
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        GeneralCosine::new(
            self.m,
            [0.42, 0.50, 0.08]
                .into_iter()
                .map(|n| W::from(n).unwrap())
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
    fn blackman_20() {
        // Created with
        // >>> from scipy.signal.windows import blackman
        // >>> blackman(20)
        let expected = vec![
            -1.38777878e-17,
            1.02226199e-02,
            4.50685843e-02,
            1.14390287e-01,
            2.26899356e-01,
            3.82380768e-01,
            5.66665187e-01,
            7.52034438e-01,
            9.03492728e-01,
            9.88846031e-01,
            9.88846031e-01,
            9.03492728e-01,
            7.52034438e-01,
            5.66665187e-01,
            3.82380768e-01,
            2.26899356e-01,
            1.14390287e-01,
            4.50685843e-02,
            1.02226199e-02,
            -1.38777878e-17,
        ];
        assert_vec_eq(expected, Blackman::new(20, true).get_window());
    }

    #[track_caller]
    fn assert_vec_eq(a: Vec<f64>, b: Vec<f64>) {
        for (a, b) in a.into_iter().zip(b) {
            assert_abs_diff_eq!(a, b, epsilon = 1e-6);
        }
    }
}
