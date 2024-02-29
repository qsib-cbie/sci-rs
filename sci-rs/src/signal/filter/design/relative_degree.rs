use nalgebra::Complex;

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

#[cfg(feature = "alloc")]
pub fn relative_degree_dyn<F>(zeros: &[Complex<F>], poles: &[Complex<F>]) -> usize {
    let degree = poles.len() as isize - zeros.len() as isize;
    if degree < 0 {
        panic!("Improper transfer function. Must have at least as many poles as zeros.");
    }
    degree as usize
}
