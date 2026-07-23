# focus


- [Features](#features)
- [Installation](#installation)
- [Quick Start: Offline vs Online
  Usage](#quick-start-offline-vs-online-usage)
  - [Offline Mode (`focus_offline`)](#offline-mode-focus_offline)
  - [Online Mode (Sequential Updates)](#online-mode-sequential-updates)
  - [Offline Mode Interface](#offline-mode-interface)
  - [Online/Sequential Interface](#onlinesequential-interface)
  - [Inspection Functions](#inspection-functions)
- [Notes](#notes)
- [Usage Examples](#usage-examples)
  - [Gaussian Univariate Detection](#gaussian-univariate-detection)
  - [One-sided Detection](#one-sided-detection)
  - [Gaussian Multivariate Detection](#gaussian-multivariate-detection)
  - [Anomaly Detection and Intensity
    Filtering](#anomaly-detection-and-intensity-filtering)
  - [Exponential family models](#exponential-family-models)
  - [Flexibility: Statistics Independent of Detector
    Type](#flexibility-statistics-independent-of-detector-type)
  - [Non-parametric changepoint
    detection](#non-parametric-changepoint-detection)
  - [AutoRegressive Process (ARP) changepoint
    detection](#autoregressive-process-arp-changepoint-detection)
- [Performance Comparison: Offline vs
  Online](#performance-comparison-offline-vs-online)
- [C++ Integration](#c-integration)
- [References](#references)
- [Authors and Contributors](#authors-and-contributors)
- [License](#license)

**focus** is a high-performance R package for online changepoint
detection in univariate and multivariate data streams. The package
provides efficient C++ implementations of the **focus** and **md-focus**
algorithms with R interfaces for real-time monitoring and offline
analysis.

## Features

- **Multiple distributions**: supports a range of models for *Gaussian
  change-in-mean*, *Poisson change-in-rate*, *Gamma/Exponential
  change-in-scale*, *Bernoulli change-in-probability*, as well as
  Non-parametric detectors. New models and cost functions are easy to
  implement!
- **Univariate and multivariate**: detect changes in univariate or
  multivariate sequences
- **Known or unknown pre-change parameters**: flexible modeling of both
  the generalised likelihood-ratio test (unknown pre-change) and the
  Page–CUSUM test (known pre-change)
- **One-sided and two-sided detection**: detects increases or decreases
  in the parameters, or both
- **C++ backend**: Language agnostic backend optimized for speed and
  scalability

## Installation

You can install the stable version of **focus** from CRAN with:

``` r
install.packages("focus")
```

You can install the development version of **focus** from source with:

``` r
# If you have devtools installed:
devtools::install_github("gtromano/unified_focus", subdir = "focus")
```

Or, if you have the package source directory:

``` r
install.packages("path/to/focus", repos = NULL, type = "source")
```

## Quick Start: Offline vs Online Usage

The package provides two modes of operation with identical statistical
results but different performance characteristics:

### Offline Mode (`focus_offline`)

All cycles and updates are handled internally in C++ for maximum
efficiency. This approach is ideal for:

- Benchmarking and performance testing
- Batch processing of complete datasets
- Computing full statistic trajectories

**Key behaviors:**

- By default, stops immediately when threshold is exceeded
- Use `threshold = Inf` to compute statistics for all observations
  (useful for visualization)

``` r
# Generate data with a changepoint
set.seed(123)
Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 2))

# Run offline detection (all processing in C++)
res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

# Plot results
plot(res$stat, type = "l", main = "FOCuS Detection Statistic (Offline)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = res$threshold, col = "red", lty = 2, lwd = 2)
if (!is.null(res$detection_time)) {
  abline(v = res$detection_time, col = "blue", lty = 2, lwd = 2)
  legend("topleft", 
         legend = c("Threshold", "Detection time"),
         col = c("red", "blue"), lty = 2, lwd = 2)
}
```

![](generate_README_files/figure-commonmark/unnamed-chunk-2-1.png)

### Online Mode (Sequential Updates)

An **online implementation** is also available—allowing you to update
the detector sequentially from R. This will be inherently slower, since
each update involves a call from R into C++. The online interface is
useful for:

- Real-time or streaming scenarios
- Integration with other R workflows
- Custom stopping rules (adaptive thresholds)
- Implementation of your own costs function

``` r
# Create detector
det <- detector_create(type = "univariate")

# Update sequentially
stat_trace <- numeric(length(Y))
threshold <- 20

for (i in seq_along(Y)) {
  detector_update(det, Y[i])
  result <- get_statistics(det, family = "gaussian")
  
  # note that the online interface supports modern R piping
  # result <- det |> detector_update(Y[i]) |> get_statistics(family = "gaussian")
  
  stat_trace[i] <- result$stat
  
  if (result$stat > threshold) {
    cat("Detection at time", i, "with changepoint estimate τ =", result$changepoint, "\n")
    stat_trace <- stat_trace[1:i]  # Truncate
    break
  }
}
```

    Detection at time 510 with changepoint estimate τ = 500 

``` r
# Plot results
plot(stat_trace, type = "l", main = "FOCuS Detection Statistic (Online)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = threshold, col = "red", lty = 2, lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-3-1.png)

Both approaches produce the same statistical results. See the
performance comparison section below for runtime benchmarks. \##
Available Functions

### Offline Mode Interface

- **`focus_offline(Y, threshold, type, family, ...)`**  
  Run the FOCuS detector in batch/offline mode with all cycles handled
  in C++ for maximum efficiency. Stops at detection by default; use
  `threshold = Inf` to compute statistics for all observations.
  - `Y`: Observation data (vector or matrix)
  - `threshold`: Detection threshold(s). Can be:
    - Scalar: Single threshold applied to all statistics
    - Vector: (In case of multiple values returned per statistics, see
      Notes below)
  - `type`: One of `"univariate"`, `"univariate_one_sided"`,
    `"multivariate"`, or `"npfocus"`
  - `family`: Distribution family - `"gaussian"`, `"poisson"`,
    `"bernoulli"`, `"gamma"`, or `"npfocus"`
  - `theta0`: (Optional) Baseline parameter vector for cost computation
  - `shape`: (Optional) Shape parameter for `family = "gamma"` (required
    for gamma)
  - `dim_indexes`: (Optional) List of dimension index vectors for
    multivariate projections
  - `quantiles`: (Optional) Quantile vector for `type = "npfocus"`
  - `pruning_mult`, `pruning_offset`: Pruning parameters (default: 2, 1)
  - `side`: Pruning side - `"right"` or `"left"` (default: `"right"`)
  - Returns: List with `stat` (matrix where each column is a statistic),
    `changepoint`, `detection_time`, `detected_changepoint`,
    `candidates`, `threshold`, `n`, `type`, `family`, and `shape` (if
    gamma)

### Online/Sequential Interface

- **`detector_create(type, ...)`**  
  Create a new online detector object. Returns an Info object pointer.
  - `type`: One of `"univariate"`, `"univariate_one_sided"`,
    `"multivariate"`, or `"npfocus"`
  - `dim_indexes`: (Optional) List of dimension index vectors for
    multivariate projections
  - `quantiles`: (Optional) Quantile vector for `type = "npfocus"`
  - `pruning_mult`, `pruning_offset`: Pruning parameters (default: 2, 1)
  - `side`: Pruning side - `"right"` or `"left"` (default: `"right"`)
- **`detector_update(detector, y)`**  
  Update the detector with a new observation vector `y`.
  - `detector`: Info object pointer from `detector_create()`
  - `y`: Observation vector (length must match detector dimension)
- **`get_statistics(detector, family, theta0 = NULL, shape = NULL)`**  
  Compute changepoint statistics for the current state.
  - `detector`: Info object pointer
  - `family`: Distribution family - `"gaussian"`, `"poisson"`,
    `"bernoulli"`, `"gamma"`, or `"npfocus"`
  - `theta0`: (Optional) Baseline parameter vector for cost computation
  - `shape`: (Optional) Shape parameter for `family = "gamma"` (required
    for gamma)
  - Returns: List with `stopping_time`, `changepoint`, and `stat`
    - `stat` can be a scalar (for single statistic families) or vector
      (for multi-statistic families like npfocus)

### Inspection Functions

- **`detector_info_n(detector)`** - Get number of observations processed
- **`detector_info_sn(detector)`** - Get cumulative sum state (vector
  for multivariate)
- **`detector_cands_len(detector)`** - Get number of candidate
  changepoints
- **`detector_candidates(detector)`** - Get all candidate changepoints
  as a data frame
  - Returns: Data frame with columns for candidate positions and their
    sufficient statistics

## Notes

- **Multiple Statistics**: Some detectors (e.g., `family = "npfocus"`)
  return multiple statistics. In `focus_offline()`, the `stat` return
  value is a matrix where each row corresponds to a time point and each
  column corresponds to a statistic.
  - For single-statistic families, the matrix has one column (R treats
    this as a vector in many contexts)
  - For multi-statistic families, use vectorized thresholds or a single
    threshold (with warning)
- **Gamma Family**: When using `family = "gamma"`, you must provide a
  positive `shape` parameter. The gamma cost function assumes this shape
  is known.

## Usage Examples

### Gaussian Univariate Detection

#### Pre-change Parameter Unknown

When the pre-change mean is unknown, FOCuS estimates it from the data:

``` r
# Generate data: 1000 obs at mean 0, then 500 obs at mean -1
set.seed(45)
Y <- c(rnorm(1000, mean = 0), rnorm(500, mean = -1))

# Offline detection (stops at detection)
res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

cat("Detection time:", res$detection_time, "\n")
```

    Detection time: 1035 

``` r
cat("Estimated changepoint:", res$detected_changepoint, "\n")
```

    Estimated changepoint: 991 

``` r
# Plot
plot(res$stat, type = "l", main = "FOCuS: Unknown Pre-change Mean",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = res$threshold, col = "red", lty = 2, lwd = 2)
abline(v = res$detection_time, col = "blue", lty = 2, lwd = 2)
abline(v = 1000, col = "green", lty = 3, lwd = 2)
legend("topleft", 
       legend = c("Statistic", "Threshold", "Detection", "True changepoint"),
       col = c("black", "red", "blue", "green"), 
       lty = c(1, 2, 2, 3), lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-4-1.png)

#### Pre-change Parameter Known

When you know the pre-change mean, provide it via `theta0` (this will
give more power to detect changes more quickly!):

``` r
set.seed(45)
theta0 <- 0
Y <- c(rnorm(1000, mean = theta0), rnorm(500, mean = theta0 - 1))

# With known pre-change parameter
res_known <- focus_offline(Y, threshold = Inf, type = "univariate", 
                           family = "gaussian", theta0 = theta0)

# Compare with unknown case
res_unknown <- focus_offline(Y, threshold = Inf, type = "univariate", 
                             family = "gaussian")

# Plot comparison
par(mfrow = c(1, 2))
plot(res_known$stat, type = "l", main = "Known θ₀ = 0",
     xlab = "Time", ylab = "Statistic", lwd = 2, ylim = range(c(res_known$stat, res_unknown$stat)))
abline(v = 1000, col = "green", lty = 3, lwd = 2)

plot(res_unknown$stat, type = "l", main = "Unknown θ₀",
     xlab = "Time", ylab = "Statistic", lwd = 2, ylim = range(c(res_known$stat, res_unknown$stat)))
abline(v = 1000, col = "green", lty = 3, lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-5-1.png)

``` r
par(mfrow = c(1, 1))
```

### One-sided Detection

Detect only increases (right-sided) or decreases (left-sided):

``` r
set.seed(789)
Y_increase <- c(rnorm(800, mean = 0), rnorm(400, mean = 1.5))

# Right-sided: detect only increases
res_right <- focus_offline(Y_increase, threshold = 30, 
                           type = "univariate_one_sided", 
                           family = "gaussian", side = "right")

cat("Right-sided detection at:", res_right$detection_time, "\n")
```

    Right-sided detection at: 816 

``` r
# Left-sided: detect only decreases (won't detect in this data)
res_left <- focus_offline(Y_increase, threshold = 30, 
                          type = "univariate_one_sided", 
                          family = "gaussian", side = "left")

cat("Left-sided detection:", 
    ifelse(is.null(res_left$detection_time), "None", res_left$detection_time), "\n")
```

    Left-sided detection: None 

``` r
# Plot comparison
par(mfrow = c(1, 2))
plot(res_right$stat, type = "l", main = "Right-sided (stat increases)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = 30, col = "red", lty = 2)
if (!is.null(res_right$detection_time)) {
  abline(v = res_right$detection_time, col = "blue", lty = 2)
}

plot(res_left$stat, type = "l", main = "Left-sided (no detection)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = 30, col = "red", lty = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-6-1.png)

``` r
par(mfrow = c(1, 1))
```

### Gaussian Multivariate Detection

For vector-valued observations:

``` r
set.seed(42)
n <- 1500
p <- 3

# Create data: changepoint at t=1000
Y_multi <- rbind(
  matrix(rnorm(1000 * p, mean = 0), ncol = p),
  matrix(rnorm(500 * p, mean = 1.2), ncol = p)
)

# Run detection
res_multi <- focus_offline(Y_multi, threshold = 30, 
                           type = "multivariate", family = "gaussian")

cat("Detection time:", res_multi$detection_time, "\n")
```

    Detection time: 1012 

``` r
cat("Estimated changepoint:", res_multi$detected_changepoint, "\n")
```

    Estimated changepoint: 1000 

``` r
cat("True changepoint: 1000\n")
```

    True changepoint: 1000

``` r
# Plot
plot(res_multi$stat, type = "l", main = "Multivariate FOCuS (p=3)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = res_multi$threshold, col = "red", lty = 2, lwd = 2)
abline(v = res_multi$detection_time, col = "blue", lty = 2, lwd = 2)
abline(v = 1000, col = "green", lty = 3, lwd = 2)
legend("topleft", 
       legend = c("Threshold", "Detection", "True changepoint"),
       col = c("red", "blue", "green"), lty = c(2, 2, 3), lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-7-1.png)

**Note:** When you have more than 5 dimensions, computing the convex
hull can become quite slow and you may not prune many points. A
practical approach is to use a low-dimensional approximation of the
hull, which often gives very similar results at a fraction of the
computation time.

Here’s an example comparing full vs low-dimensional hull computation:

``` r
set.seed(42)
n <- 1000
p <- 6

# Create data: changepoint at t=1000
Y_multi <- rbind(
  matrix(rnorm(5000 * p, mean = -1, 1), ncol = p),
  matrix(rnorm(500 * p, mean = 1.2), ncol = p)
)

# Full multivariate detection
system.time(
  res_multi <- focus_offline(Y_multi, threshold = Inf,
                             type = "multivariate", family = "gaussian")
)
```

       user  system elapsed 
      8.425   0.111   8.538 

``` r
# Low-dimensional projection approximation
dim_indexes <- generate_projection_indexes(6, 2)
system.time(
  res_multi_approx <- focus_offline(Y_multi, threshold = Inf,
                                    type = "multivariate", family = "gaussian",
                                    dim_indexes = dim_indexes)
)
```

       user  system elapsed 
      0.157   0.000   0.157 

``` r
# Verify similarity
all.equal(res_multi$stat, res_multi_approx$stat)
```

    [1] "Mean relative difference: 0.003782673"

### Anomaly Detection and Intensity Filtering

The `anomaly_intensity` parameter focuses on detecting anomalies
(epidemic changepoints) and more transient changes, while ignoring
longer, less intense changes. This is particularly useful when we have
an expected background rate, and we seek any deviations from it.

When `anomaly_intensity` is set to a positive value, candidates are
retained only if they show sufficient “signal intensity” — i.e., the
magnitude of the change relative to the segment length is at least
`anomaly_intensity`. This filtering occurs during candidate pruning and
helps reduce the number of spurious changepoints.

The following example runs the detector sequentially, and flags an alarm
whilist we’re in the anomalous period.

``` r
# Create synthetic data with brief anomalies
set.seed(999)
n <- 1000
Y_anom <- c(
  rnorm(n/2, mean = 0),
  rnorm(10, mean = -3),      # Brief, weak anomaly
  rnorm(n/2, mean = 0),
  rnorm(10, mean = 5),      # Stronger brief anomaly
  rnorm(n/2, mean = 0)
)

# Online (sequential) detection
det <- detector_create(type = "univariate", anomaly_intensity = 2)
threshold <- 20
in_anom <- FALSE
starts <- integer(0)
ends <- integer(0)

stat_trace <- numeric(length(Y_anom))

for (i in seq_along(Y_anom)) {
  detector_update(det, Y_anom[i])
  res <- get_statistics(det, family = "gaussian", theta0 = 0)
  stat <- res$stat
  stat_trace[i] <- stat

  # robustly handle NULL/NA
  if (is.null(stat) || length(stat) == 0) stat <- 0

  if (!in_anom && stat > threshold) {
    in_anom <- TRUE
    starts <- c(starts, res$changepoint)
    cat("anomaly starting at", i, "\n")
  }

  if (in_anom && stat <= threshold) {
    in_anom <- FALSE
    ends <- c(ends, res$changepoint)
    cat("anomaly ending at", i, "\n")
  }
}
```

    anomaly starting at 503 
    anomaly ending at 513 
    anomaly starting at 1011 
    anomaly ending at 1037 

``` r
# If we were still in an anomaly at the end, close it
if (in_anom) {
  ends <- c(ends, length(Y_anom))
  cat("anomaly ending at", length(Y_anom), "(end of series)", "\n")
}

# Plot the data and mark starts/ends
par(mfrow = c(2, 1))
plot(Y_anom, type = "l", main = "Data with Detected Anomalies",
     xlab = "Time", ylab = "Value", lwd = 1.2)
if (length(starts) > 0) abline(v = starts, col = "green", lty = 2)
if (length(ends) > 0) abline(v = ends, col = "red", lty = 2)
legend("topright", legend = c("an. start", "an. end"), col = c("green", "red"), lty = 2)
plot(stat_trace, type = "l", main = "Statistics Trace",
     xlab = "Time", ylab = "Statistic", lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-9-1.png)

``` r
par(mfrow = c(2, 1))
```

- **Default behavior** (`anomaly_intensity = NULL`): All candidates are
  retained based on standard pruning rules.
- **When set to a positive value**: Candidates are filtered based on the
  infinity norm of the signal-to-length ratio. A candidate survives only
  if at least one dimension has a signal magnitude ≥
  `anomaly_intensity`.

### Exponential family models

This section collects examples for the exponential-family models
implemented in `focus`: **Poisson**, **Bernoulli**, and **Gamma**. All
three examples show how to run an offline detection using
`focus::focus_offline()` – but the online, sequential iplementation
would work as well (see below for an example).

#### Bernoulli (change in success probability)

Univariate and multivariate Bernoulli examples. Here each observation is
a 0/1 draw and the success probability changes at the halfway point.

``` r
# univariate bernoulli example
set.seed(123)
n <- 2000
Y_bern <- c(rbinom(n/2, size = 1, prob = 0.2), rbinom(n/2, size = 1, prob = 0.5))

system.time({
  res_bern <- focus::focus_offline(Y_bern, threshold = Inf, type = "univariate", family = "bernoulli")
})
```

       user  system elapsed 
      0.002   0.000   0.003 

``` r
plot(res_bern$stat, main = "Bernoulli (univariate): change in success probability")
```

![](generate_README_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
# multivariate bernoulli example (two binary streams)
Y_bern_multi <- cbind(
  c(rbinom(n/2, size = 1, prob = 0.2), rbinom(n/2, size = 1, prob = 0.5)),
  c(rbinom(n/2, size = 1, prob = 0.3), rbinom(n/2, size = 1, prob = 0.6))
)

system.time({
  res_bern_multi <- focus::focus_offline(Y_bern_multi, threshold = Inf, type = "multivariate", family = "bernoulli")
})
```

       user  system elapsed 
      0.023   0.000   0.022 

``` r
plot(res_bern_multi$stat, main = "Bernoulli (multivariate): two streams")
```

![](generate_README_files/figure-commonmark/unnamed-chunk-10-2.png)

#### Poisson

A simple Poisson change-in-rate example: counts switch from a low to a
higher rate halfway through the sequence.

``` r
# Poisson change in rate example
set.seed(101)
n <- 2000
Y_pois <- c(rpois(n/2, lambda = 2), rpois(n/2, lambda = 6))

system.time({
  res_pois <- focus::focus_offline(Y_pois, threshold = Inf, type = "univariate", family = "poisson")
})
```

       user  system elapsed 
      0.002   0.000   0.002 

``` r
plot(res_pois$stat, main = "Poisson: change in rate (lambda)")
```

![](generate_README_files/figure-commonmark/unnamed-chunk-11-1.png)

#### Gamma (change in scale / rate)

The Gamma example requires specifying the `shape` parameter (positive
scalar). In the `focus` API the `shape` argument is passed to the cost
function; `theta0` is used as the null value for the relevant Gamma
parameter (scale/rate depending on your parameterisation in the
implementation). Make sure `shape > 0`. **Note**: the `shape` parameter
is specific to the **Gamma** family. If you pass `shape` when using
another family it will be ignored (hopefully).

``` r
# gamma change in scale example
set.seed(124)
n <- 2000
# shape = 2; first half scale=2, second half scale=0.5
Y_gamma <- c(rgamma(n/2, shape = 2, scale = 2), rgamma(n/2, shape = 2, scale = 0.5))

system.time({
  # note: shape = 2 must be provided for family='gamma'
  res_gamma <- focus::focus_offline(Y_gamma, threshold = Inf, type = "univariate", family = "gamma", shape = 2, theta0 = 2)
})
```

       user  system elapsed 
      0.003   0.000   0.003 

``` r
plot(res_gamma$stat, main = "Gamma: change in scale (shape = 2)")
```

![](generate_README_files/figure-commonmark/unnamed-chunk-12-1.png)

### Flexibility: Statistics Independent of Detector Type

A key feature of the **focus** package is that the detector type (how
candidates are managed) is completely independent of the statistical
model (how costs are computed). This means you can create a detector
with one data type and compute statistics using different distributional
assumptions.

Here’s an example using Poisson-generated data but comparing Gaussian
and Poisson statistics:

``` r
set.seed(2024)
# Generate Poisson count data with a rate change
Y_counts <- c(rpois(500, lambda = 10), rpois(500, lambda = 15))

# Create a univariate detector (same detector for both)
det <- detector_create(type = "univariate")

# Update with all data
for (i in seq_along(Y_counts)) {
  detector_update(det, Y_counts[i])
}

# Compare Gaussian vs Poisson statistics on the SAME detector state
result_gaussian <- get_statistics(det, family = "gaussian")
result_poisson <- get_statistics(det, family = "poisson", theta0 = 10)

cat("Using Gaussian statistic:\n")
```

    Using Gaussian statistic:

``` r
cat("  Changepoint:", result_gaussian$changepoint, "\n")
```

      Changepoint: 500 

``` r
cat("  Statistic:", round(result_gaussian$stat, 2), "\n\n")
```

      Statistic: 6426.22 

``` r
cat("Using Poisson statistic (more appropriate for count data):\n")
```

    Using Poisson statistic (more appropriate for count data):

``` r
cat("  Changepoint:", result_poisson$changepoint, "\n")
```

      Changepoint: 500 

``` r
cat("  Statistic:", round(result_poisson$stat, 2), "\n")
```

      Statistic: 521.28 

``` r
# Compute full trajectories for comparison
det2 <- detector_create(type = "univariate")
stat_gaussian <- numeric(length(Y_counts))
stat_poisson <- numeric(length(Y_counts))

for (i in seq_along(Y_counts)) {
  detector_update(det2, Y_counts[i])
  stat_gaussian[i] <- get_statistics(det2, family = "gaussian")$stat
  stat_poisson[i] <- get_statistics(det2, family = "poisson", theta0 = 10)$stat
}

# Plot comparison
par(mfrow = c(1, 2))
plot(stat_gaussian, type = "l", main = "Gaussian Statistic on Poisson Data",
     xlab = "Time", ylab = "Statistic", lwd = 2, col = "blue")
abline(v = 500, col = "green", lty = 3, lwd = 2)

plot(stat_poisson, type = "l", main = "Poisson Statistic on Poisson Data",
     xlab = "Time", ylab = "Statistic", lwd = 2, col = "red")
abline(v = 500, col = "green", lty = 3, lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-13-1.png)

``` r
par(mfrow = c(1, 1))
```

This flexibility allows you to:

- Test different statistical models on the same data
- Use the same detector infrastructure across different distributions
- Implement your own functions! The optimal change candidates are always
  accessible via:

``` r
head(detector_candidates(det2))
```

      tau   st  side
    1   0    0 right
    2 437 4263 right
    3 459 4481 right
    4 461 4502 right
    5 500 4916 right
    6 521 5210 right

### Non-parametric changepoint detection

Additionally, we find the NPFOCUS implementation for non-parametric
changepoint detection. This requires specifying quantiles to be used in
the cost function. Generally this quantiles can be estimated from the
data or specified based on domain knowledge.

Note that NPFOCuS returns two statistics: one based on summing over
quantiles, and one based on taking the maximum over quantiles. In the
offline version both statistics are returned in a matrix with two
columns.

``` r
set.seed(123)
Y <- c(rnorm(1000), rcauchy(200))

quants <- qnorm(seq(0.01, .99, length.out = 5))

res <- focus_offline(
  Y = Y,
  threshold = c(80, 25),         # detection threshold (one for sum, one for max)
  type = "npfocus",         # creates a NonparametricInfo
  family = "npfocus",       # uses compute_costs_npfocus
  quantiles = quants,       # REQUIRED for type == "npfocus"
)
par(mfrow = c(3, 1))
plot(Y[1:nrow(res$stat)])
plot(res$stat[, 1], type = "l", main = "NPFOCuS Detection Statistic (sum)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
plot(res$stat[, 2], type = "l", main = "NPFOCuS Detection Statistic (max)",
     xlab = "Time", ylab = "Statistic", lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-15-1.png)

``` r
par(mfrow = c(1, 1))
```

In the online setting, the quantile vector must be provided when
creating the detector. The statistic can be retrieved as usual, however
it will be a vector of length 2 (sum and max statistics).

``` r
set.seed(123)

det <- detector_create("npfocus", quantiles = quants)

det <- det |> 
  detector_update(Y[1]) |>
  detector_update(Y[2]) |>
  detector_update(Y[3])

det |> get_statistics(family="npfocus")
```

    $stopping_time
    [1] 3

    $changepoint
    NULL

    $stat
    [1] 3.819085 1.909543

### AutoRegressive Process (ARP) changepoint detection

The library also supports changepoint detection in AutoRegressive (AR)
processes. This is useful for detecting changes in the mean of AR(p)
processes while accounting for the temporal dependencies. The AR
coefficients can be specified (or estimated from data). As in the iid
cases, the pre-change mean can be provided if known, however this has to
be done when creating the detector (as the parameter is necessary for
pruning logic).

``` r
set.seed(123)

# Generate AR(2) process with changepoint
n_pre <- 500
n_post <- 500
ar_coefs <- c(0.7, -0.3)  # AR(2) coefficients

Y_pre <- arima.sim(n = n_pre, model = list(ar = ar_coefs), sd = 1)
Y_post <- 2 + arima.sim(n = n_post, model = list(ar = ar_coefs), sd = 1)
Y <- c(Y_pre, Y_post)

# Estimate AR parameters from pre-change data (in practice, use historical data)
ar_model <- ar(arima.sim(n = 300, model = list(ar = ar_coefs), sd = 1), order.max = 2, method = "mle")
rho_est <- ar_model$ar  # Estimated AR coefficients

# Run offline ARP detection
res <- focus_offline(
  Y = Y,
  threshold = 20,
  type = "arp",
  rho = rho_est,
  mu0_arp = 0
)

# Plot results
par(mfrow = c(2, 1))
plot(Y, type = "l", main = "AR(2) Process with Mean Shift",
     xlab = "Time", ylab = "Value")
abline(v = n_pre, col = "red", lty = 2, lwd = 2)  # True changepoint

plot(res$stat[, 1], type = "l", main = "ARP Detection Statistic",
     xlab = "Time", ylab = "Statistic", lwd = 2)
abline(h = res$threshold, col = "blue", lty = 2)
if (!is.null(res$detection_time)) {
  abline(v = res$detection_time, col = "red", lty = 2, lwd = 2)
}
```

![](generate_README_files/figure-commonmark/unnamed-chunk-17-1.png)

``` r
par(mfrow = c(1, 1))

# Show detection results
cat("Detection time:", res$detection_time, "\n")
```

    Detection time: 510 

``` r
cat("Estimated changepoint:", res$detected_changepoint, "\n")
```

    Estimated changepoint: 498 

``` r
cat("True changepoint:", n_pre, "\n")
```

    True changepoint: 500 

And in the online setting:

``` r
set.seed(123)

# Create online ARP detector
det <- detector_create("arp", rho = rho_est, mu0_arp = 0)

stat_trace <- numeric(length(Y))

for (i in seq_along(Y)) {
  detector_update(det, Y[i])
  result <- get_statistics(det, family = "arp")
  stat_trace[i] <- result$stat

  if (result$stat > 20) {
    cat("Detection at time", i, "with changepoint estimate τ =", result$changepoint, "\n")
    stat_trace <- stat_trace[1:i]  # trucate the trace
    break
  }    
    

}
```

    Detection at time 511 with changepoint estimate τ = 498 

``` r
# Plot results
plot(stat_trace, type = "l", main = "Online ARP Detection Statistic", 
     xlab = "Time", ylab = "Statistic", lwd = 2)
```

![](generate_README_files/figure-commonmark/unnamed-chunk-18-1.png)

## Performance Comparison: Offline vs Online

Here’s a direct runtime comparison between offline and online modes on
the same data:

``` r
# Generate larger dataset for benchmarking
set.seed(999)
n <- 100000
Y_bench <- c(rnorm(n/2, mean = 0), rnorm(n/2, mean = 1))

# Benchmark offline mode
cat("Offline mode (C++ loop):\n")
```

    Offline mode (C++ loop):

``` r
time_offline <- system.time({
  res_offline <- focus_offline(Y_bench, threshold = Inf, 
                               type = "univariate", family = "gaussian")
})
print(time_offline)
```

       user  system elapsed 
      0.161   0.001   0.161 

``` r
# Benchmark online mode
cat("\nOnline mode (R loop with C++ calls):\n")
```


    Online mode (R loop with C++ calls):

``` r
time_online <- system.time({
  det <- detector_create(type = "univariate")
  stat_online <- numeric(n)
  for (i in seq_along(Y_bench)) {
    detector_update(det, Y_bench[i])
    result <- get_statistics(det, family = "gaussian")
    stat_online[i] <- result$stat
  }
})
print(time_online)
```

       user  system elapsed 
       0.37    0.00    0.37 

``` r
# Verify both produce identical results
cat("\nResults identical:", all.equal(as.vector(res_offline$stat), stat_online), "\n")
```


    Results identical: TRUE 

``` r
# Speedup factor
speedup <- time_online["elapsed"] / time_offline["elapsed"]
cat("Offline mode is", round(speedup, 1), "x faster\n")
```

    Offline mode is 2.3 x faster

## C++ Integration

If you wish to use the library entirely in C++ (for maximum speed or
integration into other C++ projects), you can do so by following the
patterns in the source code. The key classes are:

- `Info` and derived classes (`UnivariateInfo`, `MultivariateInfo`)
- Cost functions (`compute_costs_gaussian`, `compute_costs_poisson`)
- `ChangepointResult` structure

Example C++ usage:

``` cpp
#include "Info.h"
#include "Costs.h"

// Create detector
auto info = std::make_shared<UnivariateInfo>();

// Update with data
for (const auto& y : data) {
    info->update({y});
    auto result = compute_costs_gaussian(*info, {theta0});
    if (result.stat.value() > threshold) {
        // Detection!
        break;
    }
}
```

## References

The software paper is available on arXiv:
<https://arxiv.org/abs/2607.19961>.

<div id="ref-pishchagina2023online" class="csl-entry">

Romano, Gaetano, Kes Ward, Yuntang Fan, Guillem Rigaill, Vincent Runge,
Idris A. Eckley, and Paul Fearnhead. 2026. “focus and focus-cpt: Fast
Online Changepoint Detection in R and Python.” *arXiv preprint
arXiv:2607.19961*. <https://arxiv.org/abs/2607.19961>

</div>

To reference the software in publications, please cite the following:

``` bibtex
@Article{,
  title = {focus and focus-cpt: Fast Online Changepoint Detection in R and Python},
  author = {Gaetano Romano and Kes Ward and Yuntang Fan and Guillem Rigaill and Vincent Runge and Idris A. Eckley and Paul Fearnhead},
  journal = {arXiv preprint arXiv:2602.04322},
  year = {2026},
  url = {https://arxiv.org/abs/2607.19961},
}
```

Other relevant references include:

<div id="ref-pishchagina2023online" class="csl-entry">

Pishchagina, Liudmila, Gaetano Romano, Paul Fearnhead, Vincent Runge,
and Guillem Rigaill. 2025. “Online Multivariate Changepoint Detection:
Leveraging Links with Computational Geometry.” *Journal of the Royal
Statistical Society Series B: Statistical Methodology*: qkaf046.
<https://doi.org/10.1093/jrsssb/qkaf046>

</div>

<div id="ref-romano2024" class="csl-entry">

Romano, Gaetano, Idris A. Eckley, and Paul Fearnhead. 2024. “A
Log-Linear Nonparametric Online Changepoint Detection Algorithm Based on
Functional Pruning.” *IEEE Transactions on Signal Processing* 72:
594–606. <https://doi.org/10.1109/tsp.2023.3343550>.

</div>

<div id="ref-romano2023fast" class="csl-entry">

Romano, Gaetano, Idris A Eckley, Paul Fearnhead, and Guillem Rigaill.
2023. “Fast Online Changepoint Detection via Functional Pruning CUSUM
Statistics.” *Journal of Machine Learning Research* 24 (81): 1–36.
<https://www.jmlr.org/papers/v24/21-1230.html>.

</div>

<div id="ref-ward2024constant" class="csl-entry">

Ward, Kes, Gaetano Romano, Idris Eckley, and Paul Fearnhead. 2024. “A
Constant-Per-Iteration Likelihood Ratio Test for Online Changepoint
Detection for Exponential Family Models.” *Statistics and Computing* 34
(3): 1–11.

</div>

</div>

## Authors and Contributors

- Gaetano Romano: [email](mailto:g.romano@lancaster.ac.uk) (**Author**)
  (**Maintainer**) (**Creator**) (**Translator**)

- Kes Ward: [email](mailto:k.ward4@lancaster.ac.uk) (**Author**)

- Yuntang Fan: [email](mailto:y.yuntang@lancaster.ac.uk) (**Author**)

- Guillem Rigaill: [email](mailto:guillem.rigaill@inrae.fr) (**Author**)

- Vincent Runge: [email](mailto:vincent.runge@univ-evry.fr) (**Author**)

- Idris A. Eckley: [email](mailto:i.eckley@lancaster.ac.uk) (**Author**)

- Paul Fearnhead: [email](mailto:p.fearnhead@lancaster.ac.uk)
  (**Author**)

## License

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see https://www.gnu.org/licenses/gpl-3.0.txt.

The package bundles code from code from Qhull (http://www.qhull.org/),
from C.B. Barber and The Geometry Center. Qhull is free software and may
be obtained via http from www.qhull.org. For details, see
[inst/COPYRIGHTS/qhull](inst/COPYRIGHTS/qhull)
