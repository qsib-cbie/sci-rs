///
/// GH (alpha-beta), GHK (alpha-beta-gamma) filters
///
pub use kalmanfilt::gh as gh_filter;

///
/// Multivariate or Univariant Kalman filters
///
pub use kalmanfilt::kalman::kalman_filter;

///
/// 1D Gaussian filters
///
pub use gaussfilt as gaussian_filter;

/// Digital IIR/FIR filter design  
/// Functions located in the [`Filter design` section of
/// `scipy.signal`](https://docs.scipy.org/doc/scipy/reference/signal.html#filter-design).
pub mod design;

mod ext;
mod sosfilt;

pub use ext::*;
pub use sosfilt::*;

#[cfg(feature = "alloc")]
mod lfilter_zi;
#[cfg(feature = "alloc")]
mod savgol_filter;
#[cfg(feature = "alloc")]
mod sosfilt_zi;
#[cfg(feature = "alloc")]
mod sosfiltfilt;

#[cfg(feature = "alloc")]
pub use lfilter_zi::*;
#[cfg(feature = "alloc")]
pub use savgol_filter::*;
#[cfg(feature = "alloc")]
pub use sosfilt_zi::*;
#[cfg(feature = "alloc")]
pub use sosfiltfilt::*;
