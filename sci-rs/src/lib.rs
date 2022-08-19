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
//! that do not require runtime allocation.`use_std` is a default feature
//! that is required to compile in some algorithms.
//!

#![cfg_attr(not(feature = "use_std"), no_std)]
#[cfg(any(not(feature = "use_std"), test))]
#[macro_use]
extern crate std;

pub mod signal;
pub mod stats;
