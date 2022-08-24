# sci-rs

Pure Rust translation of scipy for reliable cross-platform and cross-tier behavior. See the sibling repo, https://github.com/qsib-cbie/sciprs, for pythonic interaction with this library. `sci-rs` will prefer idiomatic Rust with feature parity to the scipy interface when possible.

Memory usage is a priority. While `use_std` is a default feature, the library will prefer implementations that do not require runtime allocations, unless noted otherwise.

### Python :: Replace scipy with sciprs

Here is an example of python that uses sosfiltfilt to bandpass a numpy array

```python3
import numpy as np
# from scipy.signal import butter, sosfiltfilt
from sciprs.signal import butter, sosfiltfilt

raw = np.array([0.0, 0.09414586007215595, 0.18745540640340155, ...])
sos = butter(4, [10, 50], btype='bandpass', output='sos', fs=1666)
np_ndarray_filtered_by_rust_code = sosfiltfilt(buttersos, raw, 0)
```

### Rust :: Call sci-rs directly


```rust
let y = Array1::from_vec(vec![0.0, 0.09414586007215595, 0.18745540640340155, ...]);
let sos_arr = Butter::new(4, FilterType::Bandpass(10, 50), FilterOutput::Sos, 1666.);
let sos: [Sos<f64>; N] = Sos::from_scipy(sos_arr);
let z: Vec<f64> = sosfiltfilt_dyn(y.iter(), &sos);
```

### Rust :: Conventions

Functions that end with `.*_st` do not require allocation. Functions that end with `.*_dyn` require allocation. Static implementations are preferred for maximizing the `no_std` feature set, but some algorithms require allocation. An example would be `sosfilt_st` and `sosfilt_dyn`, where benchmarks show negligible performance changes between static and dynamic implementations and 8.6x performance improvement over `scipy.signal.sosfilt` with the same data and parameters.

We use nightly with features `generic_const_exprs` and `generic_args_infer` to allow static functions to compute intermediate and output sizing when possible. As these features are not complete or stabilized, the `.*_dyn` implementations can be considered more stable than `.*_st`.
#### Rust Goal :: Copy paste python with macros

In order to minimize friction from using scipy python code directly in rust, sci-rs will include a macro system for python-compatible syntax to expand the sci-rs side logic that is currently interacting with pyo3 and the local python runtime. The following doesn't exist but provides a vision for the future compatibility.

```rust
let filtered_data: nalgebra::DVector<f64> = sciprs! {
    raw = np.array([0.0, 0.09414586007215595, 0.18745540640340155, 0.27909975437050305, 0.3682648115914595])
    sos = butter(4, [10, 50], btype='bandpass', output='sos', fs=1666)
    sosfiltfilt(sos, raw, 0)
};
```

## Build 

Install the nightly toolchain, `rustup toolchain install nightly`, `rust-toolchain.toml` asks for configures the toolchain choice.

In general, `cargo build`, `cargo test`, and `cargo doc` work like usual. `sci-rs` does not bind against your local python installation like `sciprs` will. When building ndarray or nalgebra, you should not need LAPACK or blas, but eventually it will be an option to accelerate usage on systems that have the option.

### System Linear Algebra Dependency

I haven't decided if we will be primarily using ndarray or nalgebra, for now both are supported. Either going back to python can be converted to a numpy ndarray with the python runtime. Between the two crates, one of the primary reasons for creating `sci-rs` is to support targets that may not have reasonable access to LAPACK or blas. The choice of library will be determined on feature compatibility rather than performance. That said, there will be many instances where the performance is great, and there will be benchmarks.

ndarray and nalgebra both have features related to matrixmultiply and lapack.
## Test

For correctness, results from scipy are directly compared with sciprs output. At the moment, this requires manual scripting and testing. A local python virtualenv is recommended for development of this crate.

### Gnuplot

Tests may be used to explore, so `gnuplot` is a required dev dependency. You can find more info at http://www.gnuplot.info

* macOS
    * `brew install gnuplot`
* linux
    * `apt install gnuplot`


Use `println!("x = {:?}", x);` to print out floating point data in a format that is copy-pastable into python:

```python
sin_wave = [0.0, 0.5877852522924731, 0.9510565162951535, 0.9510565162951536, 0.5877852522924732, 1.2246467991473532e-16, -0.587785252292473, -0.9510565162951535, -0.9510565162951536, -0.5877852522924734]
```

Use `gnuplot` to plot your data directly in the terminal or in a separate window
```rust
let mut fig = Figure::new();
fig.axes2d().lines(
    correlations
        .iter()
        .enumerate()
        .map(|(i, _)| i)
        .collect_vec(),
    correlations,
    &[Caption("Sin Wave")],
);
// Comment out the dumb terminal setting to make a test interactive!
fig.set_pre_commands("set term dumb 100 40");
fig.show().unwrap();
```


```
    1 +-----------------------------------------------------------------------------------+  
        |                 ************                                                      |
        |               **            **                                   Sin Wave ******* |
    0.8 |+           ***                **                                                 +|
        |          **                     **                                                |
    0.6 |+       **                         **                                             +|
        |       *                             *                                             |
        |      *                               **                                           |
    0.4 |+    *                                  *                                         +|
        |   **                                    *                                         |
    0.2 |+ *                                       *                                       +|
        | *                                         **                                      |
        |*                                            *                                     |
    0 |+                                             *                                   +|  
        |                                               *                                   |
        |                                                *                                  |
-0.2 |+                                                *                                +|   
        |                                                  **                               |
-0.4 |+                                                   *                             +|   
        |                                                     *                             |
        |                                                      *                            |
-0.6 |+                                                      **                         *|   
        |                                                         **                     ** |
-0.8 |+                                                          **                 **  +|   
        |                                                             **             **     |
        |                                                               *************       |
    -1 +-----------------------------------------------------------------------------------+ 
        0        1         2        3        4         5        6        7         8        9
```

### Maturin + Scipy + Sciprs + Numpy + Matplotlib

Once a method is implemented locally, your local `sciprs` crate can point to its dependency at your local `sci-rs` directory `sci-rs = { path = "../sci-rs" }`. You can `maturin develop` to build and install changes made in Rust and script comparisons directly through python.

```
(sciprs-env) jwtrueb@macbook-pro sciprs % maturin develop --release
🔗 Found pyo3 bindings
🐍 Found CPython 3.9 at /Users/jwtrueb/Desktop/workspace/sciprs/sciprs-env/bin/python
   Compiling autocfg v1.1.0
   Compiling target-lexicon v0.12.4
...
   Compiling nalgebra v0.31.1
   Compiling numpy v0.16.2
   Compiling sci-rs v0.1.2 (/Users/jwtrueb/Desktop/workspace/sci-rs/sci-rs)
   Compiling nalgebra-numpy v0.3.0 (https://github.com/qsib-cbie/nalgebra-numpy#427a8135)
   Compiling sciprs v0.1.0 (/Users/jwtrueb/Desktop/workspace/sciprs)
    Finished release [optimized] target(s) in 12.10s
📦 Built wheel for CPython 3.9 to /var/folders/n_/s24_t4gj6rnf7r2kxj376df40000gn/T/.tmphMfGAv/sciprs-0.1.0-cp39-cp39-macosx_11_0_arm64.whl
```

importlib.reload doesn't work, so you need to relaunch the interpreter each time.

```python
# load the module with the new Rust code
import sciprs

# do what you want
print(sciprs.call_new_function(1,2,3))
```


