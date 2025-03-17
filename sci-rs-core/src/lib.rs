//! Core library for sci-rs.

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(feature = "alloc")]
extern crate alloc;

use core::{error, fmt};

/// Errors raised whilst running sci-rs.
#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    /// Argument parsed into function were invalid.
    #[cfg(feature = "alloc")]
    InvalidArg {
        /// The invalid arg
        arg: alloc::string::String,
        /// Explaining why arg is invalid.
        reason: alloc::string::String,
    },
    /// Argument parsed into function were invalid.
    #[cfg(not(feature = "alloc"))]
    InvalidArg,
    /// Two or more optional arguments passed into functions conflict.
    #[cfg(feature = "alloc")]
    ConfictArg {
        /// Explaining what arg is invalid.
        reason: alloc::string::String,
    },
    /// Two or more optional arguments passed into functions conflict.
    #[cfg(not(feature = "alloc"))]
    ConfictArg,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
    }
}

impl error::Error for Error {}
