# focus

**focus** is a high-performance package for **online changepoint detection** in **univariate** and **multivariate** data streams.
It provides efficient C++ implementations of the **FOCUS**, **MD-FOCUS** and **NP-FOCUS** algorithms, with R and Python interfaces for **real-time monitoring** and **offline analysis**.

---

## Features

* **Multiple distributions** — supports a range of models for *Gaussian change-in-mean*, *Poisson change-in-rate*, *Gamma/Exponential change-in-scale*, *Bernoulli change-in-probability*, as well as Non-parametric detectors. New models and cost functions are easy to implement! 
* **Univariate and multivariate detection** — detect changes in univariate or multivariate sequences
* **Known or unknown pre-change parameters** — flexible modeling of both the likelihood-ratio test (unknown pre-change) and the Page–CUSUM test (known pre-change)
* **One-sided test statistics** *(univariate only)* — detects increases or decreases in the parameters, or both
* **Easy to access C++ backend** — Language agnostic backend optimized for speed and scalability

---

## Interfaces

The same **C++ backend** is shared across both interfaces:

* **R interface** → [focus/](focus)
  → See detailed documentation in [focus/README.md](focus/README.md)

* **Python interface** → [focus_cpt/](focus_cpt)
  → See detailed documentation in [focus_cpt/README](focus_cpt/README.md)

---

## Installation

### R

```r
# from source
devtools::install_github("gtromano/unified_focus", subdir = "focus")
```

### Python

```bash
# cloning the repo
# pip install focus_cpt
# or directly from GitHub
pip install git+https://github.com/gtromano/unified_focus.git#subdirectory=focus_cpt
```

---

## Example (R)

``` r
library(focus)
set.seed(45)
Y <- c(rnorm(1000, mean = 0), rnorm(500, mean = -1))

res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

cat("Detection time:", res$detection_time, "\n")
```

    Detection time: 1035 

``` r
cat("Estimated changepoint:", res$detected_changepoint, "\n")
```

    Estimated changepoint: 991 

## Example (Python)

``` python
import numpy as np
from focus_cpt import focus_offline
Y = np.concatenate([np.random.normal(0, 1, 1000), np.random.normal(-1, 1, 500)])

res = focus_offline(Y, threshold=20, type="univariate", family="gaussian")
print("Detection time:", res["detection_time"])
print("Estimated changepoint:", res["detected_changepoint"])
```

    Detection time: 1023
    Estimated changepoint: 1008

## Authors and Contributors 

* Gaetano Romano: [email](mailto:g.romano@lancaster.ac.uk) (**Author**) (**Maintainer**) (**Creator**) (**Translator**)

* Kes Ward: [email](mailto:k.ward4@lancaster.ac.uk) (**Author**)

* Yuntang Fan: [email](mailto:y.yuntang@lancaster.ac.uk) (**Author**)

* Guillem Rigaill: [email](mailto:guillem.rigaill@inrae.fr) (**Author**)

* Vincent Runge: [email](mailto:vincent.runge@univ-evry.fr) (**Author**)

* Paul Fearnhead: [email](mailto:p.fearnhead@lancaster.ac.uk) (**Author**)

* Idris A. Eckley: [email](mailto:i.eckley@lancaster.ac.uk) (**Author**)

### External libraries

* The R package version of this software bundles code from code from Qhull (http://www.qhull.org/), from C.B. Barber and The Geometry Center.
Qhull is free software and may be obtained via http from www.qhull.org. For details, see [focus/inst/COPYRIGHTS/qhull](focus/inst/COPYRIGHTS/qhull)

* The Python package depends on installation of the Qhull library. If not found, the user will be notified with instructions to install.
