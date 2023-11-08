#![no_std]

use sci_rs::na::U1;
use sci_rs::signal::filter::kalman_filter::KalmanFilter;

pub fn foo() {
    let _k = KalmanFilter::<f64, U1, U1, U1>::default();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        foo();
    }
}
