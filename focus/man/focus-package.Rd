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

  The detector, the statistics, the offline results and the projection
  index sets are S3 objects of class \code{"focus_detector"},
  \code{"focus_statistics"}, \code{"focus_offline"} and
  \code{"focus_projections"}, respectively, with \code{print} methods and,
  for detectors and offline results, \code{summary} methods; offline results
  also have a \code{plot} method. See
  \code{\link{focus-methods}}.

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
  ## Univariate Gaussian change in mean
  set.seed(123)
  Y <- c(rnorm(500, mean = 0), rnorm(500, mean = 2))

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
  summary(det, family = "gaussian")

  ## Offline (batch) example on the same data
  # Run offline detection (C++ loop). Use threshold=Inf to compute full trajectory.
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")
  if (!is.null(res$detection_time)) {
    cat("Detection at time:", res$detection_time, "\n")
  }
  summary(res)
  plot(res, data = Y)

  ## Multivariate detection
  set.seed(123)
  d <- 3
  Y_multi <- rbind(
    matrix(rnorm(300 * d, mean = 0), ncol = d),
    matrix(rnorm(100 * d, mean = 1), ncol = d)
  )
  det_mv <- detector_create(type = "multivariate")
  for (i in seq_len(nrow(Y_multi))) {
    detector_update(det_mv, Y_multi[i, ])
    r <- get_statistics(det_mv, family = "gaussian")
    if (!is.null(r$stat) && r$stat > 30) {
      cat("Multivariate detection at", i, "estimate tau =", r$changepoint, "\n")
      break
    }
  }
  summary(det_mv)

  ## Non-parametric detection (NP-FOCuS)
  # Heavy-tailed data with a change in location
  set.seed(123)
  Y_np <- c(rt(500, df = 2), rt(200, df = 2) + 1)
  # Quantiles of the pre-change distribution
  quants <- qt(c(0.1, 0.25, 0.5, 0.75, 0.9), df = 2)
  det_np <- detector_create(type = "npfocus", quantiles = quants)
  threshold_np <- c(50, 20)  # for the sum and the max statistics
  for (i in seq_along(Y_np)) {
    detector_update(det_np, Y_np[i])
    r <- get_statistics(det_np, family = "npfocus")
    if (!is.null(r$stat) && any(r$stat > threshold_np)) {
      cat("Non-parametric detection at", i, "\n")
      break
    }
  }
  r

  ## Autoregressive detection: AR(3) with a change in mean
  set.seed(123)
  rho <- c(0.5, -0.3, 0.2)  # AR coefficients
  Y_ar <- arima.sim(list(ar = rho), n = 700) + rep(c(0, 1.5), c(500, 200))
  det_ar <- detector_create(type = "arp", rho = rho, mu0_arp = 0)
  for (i in seq_along(Y_ar)) {
    detector_update(det_ar, Y_ar[i])
    r <- get_statistics(det_ar, family = "arp")
    if (!is.null(r$stat) && r$stat > 30) {
      cat("AR(3) detection at", i, "estimate tau =", r$changepoint, "\n")
      break
    }
  }
  det_ar
}
