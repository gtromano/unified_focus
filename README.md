# focus

**focus** is a high-performance package for **online changepoint detection** in **univariate** and **multivariate** data streams.
It provides efficient C++ implementations of the **FOCUS** and **MD-FOCUS** algorithms, with R and Python interfaces for **real-time monitoring** and **offline analysis**.

---

## Features

* **Multiple distributions** — supports *Gaussian* and *Poisson* families
* **Univariate and multivariate detection** — detect changes in scalar or vector-valued sequences
* **One-sided and two-sided tests** *(univariate only)* — detect increases, decreases, or both
* **Known or unknown pre-change parameters** — flexible modeling of both the likelihood-ratio test (unknown pre-change) and the Page–CUSUM test (known pre-change)
* **High-performance C++ backend** — optimized for speed and scalability

---

## Interfaces

The same **C++ backend** is shared across both interfaces:

* **R interface** → [focus/](focus)
  → See detailed documentation in [focus/README.md](focus/README.md)

* **Python interface** To be coming. Will be under → [focus_py/](focus_py)
  → See detailed documentation in [focus_py/README](focus_py/README)

---

## Installation

### R

```r
# from source
devtools::install_github("yourusername/focus", subdir = "focus")
```

### Python (TBC!!)

```bash
pip install .
# or directly from GitHub
pip install git+https://github.com/yourusername/focus.git#subdirectory=focus_py
```

---

## Example (R)

```r
library(focus)

# Univariate Gaussian example
set.seed(1)
x <- c(rnorm(100, 0), rnorm(100, 2))
res <- focus_offline(x, threshold = 10, type = "univariate", family = "gaussian")
plot(res$stat, type = "l", main = "FOCUS statistic over time")
abline(h = res$threshold, col = "red", lty = 2)
```
