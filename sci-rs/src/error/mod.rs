use core::{error, fmt};

/// Errors raised whilst running sci-rs.
#[derive(Debug, PartialEq, Eq)]
#[cfg(feature = "alloc")]
pub enum Error {
    /// Argument parsed into function were invalid.
    InvalidArg {
        /// The invalid arg
        arg: alloc::string::String,
        /// Explaining why arg is invalid.
        reason: alloc::string::String,
    },
    /// Two or more optional arguments passed into functions conflict.
    ConfictArg {
        /// Explaining what arg is invalid.
        reason: alloc::string::String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        todo!()
    }
}

impl error::Error for Error {}
