//! Embeddable Digital Signal Processing
//!
//! `sci_rs` is a collection of pure Rust implementations of
//! useful statistical and filtering constructs.
//!
//! The goal for `sci_rs` is to provide f32 and f64 generic routines
//! for processing real data with the same code used in research
//! and deployment. We want predictable behavior on various compute tiers.
//!
//! Where analogous, interfaces will try to mirror scipy and numpy so
//! that exploratory analysis is performed by the same code that runs in
//! production.
//!
//! Where allocation is not necessary, this crate will support algorithms
//! that do not require runtime allocation. `alloc` is a default feature
//! that is required to compile in some algorithms.
//!

#![allow(unused)]
#![cfg_attr(not(feature = "std"), no_std)]
#![deny(missing_docs)]

#[cfg(feature = "alloc")]
extern crate alloc;

/// Re-export nalgebra for convenience and to avoid version conflicts
pub use nalgebra as na;

/// Linear algebra
pub mod linalg;

/// Digital signal processing
pub mod signal;

/// Statistics
pub mod stats;

/// Special math functions
pub mod special;

/// Debug plotting
#[cfg(feature = "plot")]
pub mod plot;
