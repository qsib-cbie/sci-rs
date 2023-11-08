pub use kalmanfilt::gh as gh_filter;
pub use kalmanfilt::kalman::kalman_filter;

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
