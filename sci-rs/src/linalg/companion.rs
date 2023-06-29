use ::core::{
    borrow::Borrow,
    ops::{Div, Neg},
};
use nalgebra::{allocator::Allocator, *};
use num_traits::{One, Zero};

#[cfg(feature = "use_std")]
pub fn companion_dyn<T, B, I>(itr: I, m: usize) -> OMatrix<T, Dyn, Dyn>
where
    T: Scalar + One + Zero + Div<Output = T> + Neg<Output = T> + Copy,
    B: Borrow<T>,
    I: Iterator<Item = B>,
    DefaultAllocator: Allocator<T, Dyn, Dyn>,
{
    let mut itr = itr;
    let a0: T = *itr.next().expect("Invalid data length").borrow();
    let itr = itr
        .enumerate()
        .map(|(i, ai)| ((0, i), -*(ai.borrow()) / a0))
        .chain((0..(m - 2)).map(|i| (((i + 1), i), T::one())));
    let mut m = Matrix::<
        T,
        Dyn,
        Dyn,
        <DefaultAllocator as allocator::Allocator<T, Dyn, Dyn>>::Buffer,
    >::zeros(m - 1, m - 1);
    for (i, t) in itr {
        unsafe {
            *m.get_unchecked_mut(i) = t;
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scipy_example_dyn() {
        const M: usize = 4;
        let data: [_; M] = [1, -10, 31, -30];
        let matrix: DMatrix<_> = companion_dyn(data.iter().map(|i| *i as f32), data.len());

        let expected = matrix!(
            10., -31.,  30.;
            1.,   0.,   0.;
            0.,   1.,   0.;
        );

        assert_eq!(expected, matrix);
    }
}
