use super::{extend, len_guard, truncate};
use num_traits::{real::Real, Float};

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `Triangle` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct Triangle {
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Whether the window is symmetric.
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl Triangle {
    /// Returns a Triangle struct.  
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, sym: bool) -> Self {
        Triangle { m, sym }
    }
}

#[cfg(feature = "alloc")]
impl<W> GetWindow<W> for Triangle
where
    W: Real,
{
    /// Return a window of type: Triangle.
    ///
    /// This is the not the same as the Bartlett window.
    ///
    /// # Parameters
    /// `self`: [Triangle]
    ///
    /// # Returns
    /// `w`: `vec<F>`  
    ///     The window, with the maximum value normalized to 1 (though the value 1 does not appear
    ///     if `M` is even and `sym` is True).
    ///
    /// # Example
    /// ```
    /// use sci_rs::signal::windows::{Triangle, GetWindow};
    ///
    /// let nx = 8;
    /// let tri = Triangle::new(nx, true);
    /// assert_eq!(vec![0.125, 0.375, 0.625, 0.875, 0.875, 0.625, 0.375, 0.125], tri.get_window());
    ///
    /// let nx = 9;
    /// let tri = Triangle::new(nx, false);
    /// assert_eq!(vec![0.1, 0.3, 0.5, 0.7, 0.9, 0.9, 0.7, 0.5, 0.3], tri.get_window());
    /// ```
    ///
    /// # References
    /// <https://en.wikipedia.org/wiki/Window_function#Triangular_window>
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        if len_guard(self.m) {
            return Vec::<W>::new();
        }
        let (m, needs_trunc) = extend(self.m, self.sym);

        let mut n = (1..=((m + 1) / 2)).map(|x| W::from(x).unwrap());
        let m_f: W = W::from(m).unwrap();
        let w: Vec<W> = match m % 2 {
            0 => {
                let mut w: Vec<W> = n
                    .map(|n| (W::from(2).unwrap() * n - W::one()) / m_f)
                    .collect();
                w.extend(w.clone().iter().rev());
                w
            }
            1 => {
                let mut w: Vec<W> = n
                    .map(|n| W::from(2).unwrap() * n / (m_f + W::one()))
                    .collect();
                w.extend(w.clone().iter().rev().skip(1));
                w
            }
            _ => panic!(),
        };

        truncate(w, needs_trunc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    #[cfg(feature = "alloc")]
    fn case_even_true() {
        // from scipy.signal.windows import triangle
        // triangle(n)

        let upper = 1_000;
        for i in 0..upper {
            let nx = 2 * i;
            let tri = Triangle::new(nx, true);
            let expected: Vec<f64> = (0..nx)
                .into_iter()
                .chain((0..nx).rev())
                .filter(|n| n % 2 == 1)
                .map(|n| n as f64 / nx as f64)
                .collect();
            let result: Vec<f64> = tri.get_window();

            assert_eq!(expected, result);
            // for (&e, t) in expected.iter().zip(tri.get_window::<f64>()) {
            //     assert_relative_eq!(e, t, max_relative = 1e-7)
            // }
        }
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn case_even_false() {
        // from scipy.signal import get_window
        // get_window('triangle', 4)

        let upper = 1_000;
        for i in 0..upper {
            let nx = 2 * i;
            let tri = Triangle::new(nx, false);
            let expected: Vec<f64> = (1..=(nx / 2) + 1)
                .into_iter()
                .chain((1..(nx / 2) + 1).rev())
                .take(nx)
                .map(|n| n as f64 / ((nx / 2) + 1) as f64)
                .collect();

            let result: Vec<f64> = tri.get_window();
            assert_eq!(expected, result);
            // for (&e, t) in expected.iter().zip(tri.get_window::<f64>()) {
            //     assert_relative_eq!(e, t, max_relative = 1e-7)
            // }
        }
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn case_odd_true() {
        // from scipy.signal.windows import triangle
        // triangle(nx)

        let upper = 1_000;
        for i in 1..upper {
            let nx = 2 * i + 1;
            let tri = Triangle::new(nx, true);
            let expected: Vec<_> = (1..=(nx + 1) / 2)
                .into_iter()
                .chain((1..=(nx + 1) / 2).rev().skip(1))
                .map(|n| (n as f64) / ((nx + 1) as f64 / 2.))
                .collect();

            let result: Vec<f64> = tri.get_window();
            assert_eq!(expected, result);
            // for (&e, t) in expected.iter().zip(tri.get_window::<f64>()) {
            //     assert_relative_eq!(e, t);
            // }
        }
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn case_odd_false() {
        // from scipy.signal import get_window
        // get_window('triangle', 5)

        let upper = 1_000;
        for i in 1..upper {
            let nx = 2 * i + 1;
            let tri = Triangle::new(nx, false);
            let expected: Vec<_> = (0..=nx)
                .into_iter()
                .chain((0..=nx).into_iter().rev())
                .filter(|n| n % 2 == 1)
                .take(nx)
                .map(|n| n as f64 / (nx as f64 + 1.))
                .collect();

            let result: Vec<f64> = tri.get_window();
            assert_eq!(expected, result);
            // for (&e, t) in expected.iter().zip(tri.get_window::<f64>()) {
            //     assert_relative_eq!(e, t);
            // }
        }
    }
}
