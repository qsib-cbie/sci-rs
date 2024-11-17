use std::io::Write;

/// Debug utility function that will run a python script to plot the data.
///
/// This function generates a Python script to create plots of the input data and their autocorrelations.
/// It then executes the script using the system's Python interpreter.
///
/// Note: This function will open a new window to display the plots and will block execution until the window is closed.
/// It also suppresses stdout and stderr from the Python process.
pub fn python_plot(xs: Vec<&[f32]>) {
    let script = format!(
        r#"
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import correlate

xs = {:?}
fig = plt.figure(figsize=(12, 12))
gs = gridspec.GridSpec(len(xs), 2)
for i, x in enumerate(xs):
    ax = plt.subplot(gs[i, 0])
    ax.plot(x, label = f"C{{i}}")
    ax.legend()
    ax.set_xlabel("Samples")
    ax = plt.subplot(gs[i, 1])
    autocorr = correlate(x, x, mode='full')
    normcorr = autocorr / autocorr.max()
    offsets = range(-len(x) + 1, len(x))
    ax.plot(offsets, normcorr, label = f"Autocorrelation of C{{i}}")
    ax.legend()
    ax.set_xlabel("Lag")
plt.show()
"#,
        xs
    );
    // Run the script with python
    let script = script.as_bytes();
    let mut python = match std::process::Command::new("python")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::null()) // noisy
        .stderr(std::process::Stdio::null()) // noisy
        .spawn()
    {
        Ok(p) => p,
        Err(_) => return, // Return early if python fails to start
    };

    if let Some(mut stdin) = python.stdin.take() {
        if stdin.write_all(script).is_err() {
            return; // Return early if writing fails
        }
    } else {
        return; // Return early if we can't get stdin
    }

    // Wait for the python process to finish, ignoring any errors
    let _ = python.wait(); // Ignore errors as this may be called from CI
}
