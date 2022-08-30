use heapless::Vec;
use nalgebra::Complex;

pub fn relative_degree<F, const M: usize>(
    zeros: &Vec<Complex<F>, M>,
    poles: &Vec<Complex<F>, M>,
) -> usize {
    let degree = poles.len() as isize - zeros.len() as isize;
    if degree < 0 {
        panic!("Improper transfer function. Must have at least as many poles as zeros.");
    }
    degree as usize
}
