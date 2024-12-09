use super::{extend, len_guard};
use num_traits::real::Real;

#[cfg(feature = "alloc")]
use super::GetWindow;
#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

/// Collection of arguments for window `Boxcar` for use in [GetWindow].
#[derive(Debug, Clone, PartialEq)]
pub struct Boxcar {
    /// Number of points in the output window. If zero, an empty array is returned in [GetWindow].
    pub m: usize,
    /// Whether the window is symmetric. (Has no effect for boxcar.)
    ///
    /// When true, generates a symmetric window, for use in filter design.  
    /// When false, generates a periodic window, for use in spectral analysis.
    pub sym: bool,
}

impl Boxcar {
    /// Returns a Boxcar struct.
    ///
    /// # Parameters
    /// * `m`:  
    ///     Number of points in the output window. If zero, an empty array is returned.  
    /// * `sym`:  
    ///     When true, generates a symmetric window, for use in filter design.  
    ///     When false, generates a periodic window, for use in spectral analysis.
    pub fn new(m: usize, sym: bool) -> Self {
        Boxcar { m, sym }
    }
}

#[cfg(feature = "alloc")]
impl<W> GetWindow<W> for Boxcar
where
    W: Real,
{
    /// Return a window of type: boxcar or rectangular window.
    ///
    /// Also known as a rectangular window or Dirichlet window, this is equivalent to no window at
    /// all.
    ///
    /// # Parameters
    /// `self`: [Boxcar]
    ///
    /// # Example
    /// ```
    /// use sci_rs::signal::windows::{Boxcar, GetWindow};
    ///
    /// let nx = 5;
    /// let boxcar = Boxcar::new(nx, false);
    /// let window: Vec<f64> = boxcar.get_window();
    /// assert_eq!(vec![1.; nx], window);
    /// ```
    #[cfg(feature = "alloc")]
    fn get_window(&self) -> Vec<W> {
        if len_guard(self.m) {
            return Vec::<W>::new();
        }
        let (m, needs_trunc) = extend(self.m, self.sym);

        if !needs_trunc {
            vec![W::one(); m]
        } else {
            vec![W::one(); m - 1]
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn case_a() {
        let nx = 5;
        let boxcar = Boxcar::new(nx, false);
        let window: Vec<f64> = boxcar.get_window();
        assert_eq!(vec![1.; nx], window);
    }
}
