/// Digital Filtering
pub mod filter;

/// Signal Generation
pub mod wave;

/// Convolution
#[cfg(feature = "std")]
pub mod convolve;

/// Signal Resampling
#[cfg(feature = "std")]
pub mod resample;
