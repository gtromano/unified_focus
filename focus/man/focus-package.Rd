\name{focus-package}
\alias{focus-package}
\alias{focus}
\docType{package}
\title{
  Online Changepoint Detection in Univariate and Multivariate Data Streams
}
\description{
  The \pkg{focus} package provides high-performance implementations of
  online and offline changepoint detection algorithms for univariate and
  multivariate data streams. The package exposes optimized C++ backends
  (FOCuS and md-FOCuS style algorithms) with R interfaces for both
  real-time monitoring (sequential updates) and batch/offline analysis.
}
\details{
  \pkg{focus} implements efficient changepoint detection for a range of
  statistical models and use-cases.  The package supports multiple
  distributional families (Gaussian, Poisson, Bernoulli, Gamma) and a
  non-parametric NPFOCuS variant.  It provides two primary modes:

  \describe{
    \item{Offline mode (\code{focus_offline})}{Processes all observations
      in C++ for maximum efficiency. Useful for benchmarking, for
      computing full statistic trajectories, or for batch processing.
      By default the offline call stops when a threshold is exceeded;
      use \code{threshold = Inf} to compute statistics for all
      observations.}
    \item{Online / sequential mode (\code{detector_create},
      \code{detector_update}, \code{get_statistics})}{Create an online
      detector object and update it one observation at a time from R.
      This provides a flexible streaming API and access to candidate
      segments, at the cost of more R/C calls per observation.}
  }

  Main features:
  \itemize{
    \item Multiple distributions: Gaussian, Poisson, Bernoulli, Gamma,
      and non-parametric NPFOCuS.
    \item Univariate and multivariate observations.
    \item Known or unknown pre-change parameters (generalised
      likelihood-ratio and Page–CUSUM style tests).
    \item One-sided and two-sided detection.
    \item Flexible interface: statistical cost functions are independent
      from the detector candidate-management strategy.
    \item High-performance C++ backend designed for integration with R
      and direct use from C++ projects.
  }

  Important functions:
  \itemize{
    \item \code{focus_offline(Y, threshold, type, family, ...)}:
      run complete offline detection; returns full statistic trajectories,
      detected changepoints, candidate segments, thresholds and metadata.
    \item \code{detector_create(type, ...)}: construct an online detector
      object.
    \item \code{detector_update(det, y)}: update the detector with a new
      observation.
    \item \code{get_statistics(det, family, theta0 = NULL, shape = NULL)}:
      compute statistics for the current detector state.
    \item Inspection utilities: \code{detector_info_n()}, \code{detector_info_sn()},
      \code{detector_cands_len()}, \code{detector_candidates()}.
  }

}
\author{
  Gaetano Romano [aut, cre], \email{g.romano@lancaster.ac.uk}\\
  Kes Ward [aut], \email{k.ward4@lancaster.ac.uk}\\
  Yuntang Fan [aut], \email{y.yuntang@lancaster.ac.uk}\\
  Guillem Rigaill [aut], \email{guillem.rigaill@inrae.fr}\\
  Vincent Runge [aut], \email{vincent.runge@univ-evry.fr}\\
  Paul Fearnhead [aut], \email{p.fearnhead@lancaster.ac.uk}\\
  Idris A. Eckley [aut], \email{i.eckley@lancaster.ac.uk}\\
  Maintainer: Gaetano Romano \email{g.romano@lancaster.ac.uk}
}
\references{
  Pishchagina, L., G. Romano, P. Fearnhead, V. Runge, and G. Rigaill (2025).
  Online Multivariate Changepoint Detection: Leveraging Links with
  Computational Geometry. \emph{JRSS B}. doi:10.1093/jrsssb/qkaf046.

  Romano, G., I. A. Eckley, and P. Fearnhead (2024). A Log-Linear
  Nonparametric Online Changepoint Detection Algorithm Based on
  Functional Pruning. \emph{IEEE Transactions on Signal Processing}, 72.

  Romano, G., I. A. Eckley, P. Fearnhead, and G. Rigaill (2023).
  Fast Online Changepoint Detection via Functional Pruning CUSUM
  Statistics. \emph{Journal of Machine Learning Research}, 24(81):1–36.

  Ward, K., G. Romano, I. Eckley, and P. Fearnhead (2024). A
  Constant-Per-Iteration Likelihood Ratio Test for Online Changepoint
  Detection for Exponential Family Models. \emph{Statistics and
  Computing}, 34(3).

  \emph{Note:} The package bundles and links to code from Qhull (C. B.
  Barber et al., The Geometry Center). Qhull functions are used for
  convex-hull computations in multivariate pruning. See \file{inst/COPYRIGHTS/}
  for the Qhull license and attribution details.
}
\keyword{ package }
\seealso{
  \code{\link{focus_offline}}, \code{\link{detector_create}},
  \code{\link{detector_update}}, \code{\link{get_statistics}},
  \code{\link{detector_candidates}}, \code{\link{generate_projection_indexes}}
}
\examples{
  \donttest{
    ## Offline (batch) example: univariate Gaussian change-in-mean
    set.seed(123)
    Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 2))

    # Run offline detection (C++ loop). Use threshold=Inf to compute full trajectory.
    res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")
    if (!is.null(res$detection_time)) {
      cat("Detection at time:", res$detection_time, "\n")
    }

    ## Online (sequential) example
    det <- detector_create(type = "univariate")
    stat_trace <- numeric(length(Y))
    threshold <- 20
    for (i in seq_along(Y)) {
      detector_update(det, Y[i])
      r <- get_statistics(det, family = "gaussian")
      stat_trace[i] <- r$stat
      if (!is.null(r$stat) && r$stat > threshold) {
        cat("Online detection at", i, "estimate tau =", r$changepoint, "\n")
        break
      }
    }

    ## Multivariate offline example (p = 3)
    set.seed(42)
    p <- 3
    Y_multi <- rbind(
      matrix(rnorm(1000 * p, mean = 0), ncol = p),
      matrix(rnorm(500 * p, mean = 1.2), ncol = p)
    )
    res_multi <- focus_offline(Y_multi, threshold = 30, type = "multivariate", family = "gaussian")
    cat("Multivariate detection time:", res_multi$detection_time, "\n")
  }
}
