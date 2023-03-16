#[cfg(feature = "unstable")]
use heapless::Vec;
use nalgebra::Complex;

#[cfg(feature = "unstable")]
pub fn relative_degree_st<F, const M: usize>(
    zeros: &Vec<Complex<F>, M>,
    poles: &Vec<Complex<F>, M>,
) -> usize {
    let degree = poles.len() as isize - zeros.len() as isize;
    if degree < 0 {
        panic!("Improper transfer function. Must have at least as many poles as zeros.");
    }
    degree as usize
}

#[cfg(feature = "use_std")]
pub fn relative_degree_dyn<F>(zeros: &Vec<Complex<F>>, poles: &Vec<Complex<F>>) -> usize {
    let degree = poles.len() as isize - zeros.len() as isize;
    if degree < 0 {
        panic!("Improper transfer function. Must have at least as many poles as zeros.");
    }
    degree as usize
}
