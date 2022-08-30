# Unreleased

- Added:
    - Began using more RealField and ComplexField, Float/RealField/ComplexField will change
    - `no_std` compat butterworth iirfilter design for ba, zpk, or sos
    - signal
        - design
            - `butter_st`
            - `iirfilter_st`
            - `cplxreal_st`
            - `sort_cplx_st`
            - `buttap_st`
            - `lp2bp_zpk_st`
            - `zpk2sos_st`
            - `zpk2tf_st`
            - `poly_st` (1D only)
            - `FilterOutputType`
            - `BaFormatFilter`
            - `ZpkFormatFilter`
            - `SosFormatFilter`
            - `DigitalFilter`
            - `FilterType`
            - `FilterBandType`

---

# v0.1.4

## sci-rs

- Added:
    - signal
        - design
            - `Sos`
            - `Sos::from_scipy`
        - `sosfilt_st`
            - 🔥 8.6x faster 🚀
        - `sosfilt_dyn`
            - 🔥 8.6x faster 🚀
        - `sosfilt_zi_dyn`
        - `lfilter_zi_dyn`
        - `sosfiltfilt_dyn`
            - 🔥 4.5x faster 🚀
        - `odd_ext_dyn`
    - linalg
        - `companion`
- Benchmarks
    - Generally speaking until there is a direct comparison with the Python runtime and GIL via sciprs, speedups don't directly traslate.
    - Rust
        - Criterion in Rust benchmark
    - Python
        - timeit.timeit over 100 or 1000 runs
- Dependencies changed
    - `nalgebra` with no_std optional
    - `ndarray` with no_std default
- Updated License to MIT or Apache 2.0

---

*CHANGELOG begins...*