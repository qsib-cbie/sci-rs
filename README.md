# sci-rs

A Rust implementation for scipy functionality. A pure Rust implementation for reliable cross-platform and cross-tier behavior.

See the sibling repo, https://github.com/qsib-cbie/sciprs, for pythonic interaction with this library. Otherwise, checkout the docs to use this library directly from Rust.

Memory usage is a priority. While `use_std` is a default feature, the library will prefer implementations that do not require runtime allocations, unless noted otherwise.

## Example

We sci-rs will prefer idiomatic Rust with parity to the scipy interface when possible

```rust
// Here is some floating point data
let sin_wave: Vec<f64> = todo!();

// Here is a digital filter designed in scipy
let filter: [f64; 24] = [todo!()];

// We get the same data as scipy sosfilt with this filter
let filtered_wave = sosfilt(sin_wave.iter(), &sos).collect_vec();
```

## System Build Dependencies

In general, `cargo build` and `cargo test` like usual. 

Tests may be used to explore, so `gnuplot` is a required dev dependency. You can find more info at http://www.gnuplot.info

* macOS
    * `brew install gnuplot`
* linux
    * `apt install gnuplot`

## Useful Test Interactions

Use `println!("x = {:?}", x);` to print out floating point data in a format that is copy-pastable into python:

```python
sin_wave = [0.0, 0.5877852522924731, 0.9510565162951535, 0.9510565162951536, 0.5877852522924732, 1.2246467991473532e-16, -0.587785252292473, -0.9510565162951535, -0.9510565162951536, -0.5877852522924734]
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
    &[Caption("Sin Wave Autocorrelation")],
);
// Comment out the dumb terminal setting to make a test interactive!
fig.set_pre_commands("set term dumb 120 40");
fig.show().unwrap();
```

```python
correlations = [1.0, 0.809017, 0.37811527, -0.08541019, -0.4045085, -0.5, -0.4045085, -0.2236068, -0.0690983, -7.006946e-10]
```

```                                                                                                                
                                                                                                    
       1 +-----------------------------------------------------------------------------------+      
         | ***                                                                               |      
         |    ***                                           Sin Wave Autocorrelation ******* |      
     0.8 |+      **                                                                         +|      
         |         **                                                                        |      
         |           *                                                                       |      
     0.6 |+           *                                                                     +|      
         |             **                                                                    |      
         |               *                                                                   |      
         |                **                                                                 |      
     0.4 |+                 *                                                               +|      
         |                   *                                                               |      
         |                    *                                                              |      
     0.2 |+                    *                                                            +|      
         |                      **                                                           |      
         |                        *                                                          |      
       0 |+                        *                                                     ****|      
         |                          *                                              ******    |      
         |                           *                                         ****          |      
    -0.2 |+                           **                                    ***             +|      
         |                              **                               ***                 |      
         |                                **                          ***                    |      
         |                                  **                     ***                       |      
    -0.4 |+                                   *****          ******                         +|      
         |                                         **********                                |      
         |                                                                                   |      
    -0.6 +-----------------------------------------------------------------------------------+      
         0        1         2        3        4         5        6        7         8        9      ```
