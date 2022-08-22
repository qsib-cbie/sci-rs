use ::core::{
    borrow::Borrow,
    ops::{Div, Neg},
};
use nalgebra::{allocator::Allocator, *};
use num_traits::{One, Zero};

pub fn companion_st<T, I, const M: usize>(itr: I) -> OMatrix<T, Const<{ M - 1 }>, Const<{ M - 1 }>>
where
    T: Scalar + One + Zero + Div<Output = T> + Neg<Output = T> + Copy,
    I: Iterator<Item = T>,
    DefaultAllocator: Allocator<T, Const<{ M - 1 }>, Const<{ M - 1 }>>,
{
    let mut itr = itr;
    let a0: T = *itr.next().expect("Invalid data length").borrow();
    let itr = itr
        .enumerate()
        .map(|(i, ai)| ((0, i), -*(ai.borrow()) / a0))
        .chain((0..(M - 2)).map(|i| (((i + 1), i), T::one())));
    let mut m = Matrix::<
        T,
        Const<{ M - 1 }>,
        Const<{ M - 1 }>,
        <DefaultAllocator as allocator::Allocator<T, Const<{ M - 1 }>, Const<{ M - 1 }>>>::Buffer,
    >::zeros();
    for (i, t) in itr {
        unsafe {
            *m.get_unchecked_mut(i) = t;
        }
    }
    m
}

pub fn companion_dyn<T, B, I>(itr: I, m: usize) -> OMatrix<T, Dynamic, Dynamic>
where
    T: Scalar + One + Zero + Div<Output = T> + Neg<Output = T> + Copy,
    B: Borrow<T>,
    I: Iterator<Item = B>,
    DefaultAllocator: Allocator<T, Dynamic, Dynamic>,
{
    let mut itr = itr;
    let a0: T = *itr.next().expect("Invalid data length").borrow();
    let itr = itr
        .enumerate()
        .map(|(i, ai)| ((0, i), -*(ai.borrow()) / a0))
        .chain((0..(m - 2)).map(|i| (((i + 1), i), T::one())));
    let mut m = Matrix::<
        T,
        Dynamic,
        Dynamic,
        <DefaultAllocator as allocator::Allocator<T, Dynamic, Dynamic>>::Buffer,
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
    fn scipy_example_st() {
        const M: usize = 4;
        let data: [_; M] = [1, -10, 31, -30];
        let matrix: SMatrix<_, _, _> = companion_st::<_, _, M>(data.iter().map(|i| *i as f32));

        let expected = matrix!(
            10., -31.,  30.;
            1.,   0.,   0.;
            0.,   1.,   0.;
        );

        assert_eq!(expected, matrix);
    }

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
