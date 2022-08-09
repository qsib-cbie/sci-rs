# Embedded Rust Digital Signal Processing

This library contains top level functions for processing signal. Functions are intended to handle f32 or f64 by value or reference. 

Memory usage is a priority. Without enabling compilation of std, the library is `no_std` by default and will not do any heap allocations.

## System Build Dependencies

In general, `cargo build` and `cargo test` like usual. 

Tests may be used to explore, so `gnuplot` is a required dev dependency. You can find more info at http://www.gnuplot.info


On Mac OS:
* `brew install gnuplot`

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
