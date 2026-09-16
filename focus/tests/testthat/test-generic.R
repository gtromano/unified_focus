library(testthat)
library(focus)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_gaussian <- function(n_pre = 500, n_post = 500, mu_post = 2, seed = 123) {
  set.seed(seed)
  c(rnorm(n_pre, mean = 0), rnorm(n_post, mean = mu_post))
}

focus_sequential_test <- function(
    Y,
    threshold = Inf,
    type = "univariate",
    family = "gaussian",
    theta0 = NULL,
    shape = NULL,
    side = "right",
    quantiles = NULL,
    rho = NULL,
    mu0_arp = NULL
) {

  det <- detector_create(type = type, side = side, quantiles = quantiles, rho = rho, mu0_arp = mu0_arp)

  n_obs <- if (is.matrix(Y)) nrow(Y) else length(Y)

  stat_trace <- vector("list", n_obs)
  changepoint_trace <- vector("list", n_obs)

  detection_time <- NULL
  detected_changepoint <- NULL

  for (i in seq_len(n_obs)) {
    if (is.matrix(Y)) {
      detector_update(det, Y[i, ])
    } else {
      detector_update(det, Y[i])
    }

    r <- get_statistics(det, family = family, theta0 = theta0, shape = shape)

    stat_trace[[i]] <- r$stat
    changepoint_trace[[i]] <- r$changepoint

    crossed <- FALSE
    if (!is.null(r$stat)) {
      if (is.matrix(r$stat) || (is.array(r$stat) && length(dim(r$stat)) > 1)) {
        crossed <- all(r$stat > threshold)
      } else if (length(r$stat) > 1 && length(threshold) == 1) {
        crossed <- any(r$stat > threshold)
      } else if (length(threshold) > 1) {
        crossed <- all(r$stat > threshold)
      } else {
        crossed <- isTRUE(r$stat > threshold)
      }
    }

    if (crossed) {
      detection_time <- i
      detected_changepoint <- r$changepoint
      break
    }
  }

  kept_n <- if (is.null(detection_time)) n_obs else detection_time

  stat_trace <- stat_trace[seq_len(kept_n)]
  changepoint_trace <- changepoint_trace[seq_len(kept_n)]

  stat_out <- if (kept_n == 0) {
    numeric(0)
  } else if (is.numeric(stat_trace[[1]]) && length(stat_trace[[1]]) == 1) {
    unlist(stat_trace)
  } else {
    do.call(rbind, stat_trace)
  }

  changepoint_out <- if (kept_n == 0) {
    numeric(0)
  } else {
    unlist(changepoint_trace)
  }


  list(
    stat = stat_out,
    changepoint = changepoint_out,
    detected_changepoint = detected_changepoint,
    detection_time = detection_time,
    candidates = detector_candidates(det),
    threshold = threshold,
    n = kept_n,
    type = type,
    family = family,
    shape = shape
  )
}



# ===========================================================================
# focus_sequential_test — Gaussian univariate
# ===========================================================================

test_that("sequential Gaussian: detects changepoint at correct time (seed 123)", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_sequential_test(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$detection_time,      510)
  expect_equal(res$detected_changepoint, 500)
  expect_equal(res$n,                    510)
})


test_that("sequential Gaussian: detects changepoint at correct time, shifted by 123 (seed 123)", {
  Y   <- make_gaussian(seed = 123) + 123
  res <- focus_sequential_test(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$detection_time,      510)
  expect_equal(res$detected_changepoint, 500)
  expect_equal(res$n,                    510)
})

test_that("sequential Gaussian: check the values of the statistic (theta0 = 123) (seed 123)", {
  Y   <- make_gaussian(seed = 123) + 123
  res <- focus_sequential_test(Y, threshold = 20, theta0 = 123, type = "univariate", family = "gaussian")

  expect_equal(res$stat[509],      21.04167653450661)
  expect_equal(res$stat[456],      7.080510139397964)
})

test_that("sequential Gaussian: check the values of the statistic (seed 123)", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_sequential_test(Y, threshold = 20,  type = "univariate", family = "gaussian")

  expect_equal(res$stat[510],      23.05886676141565)
  expect_equal(res$stat[509],      19.74501146734024)
  expect_equal(res$stat[456],      7.158358180852482)

})

test_that("sequential Gaussian: stat vector length equals detection time", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_sequential_test(Y, threshold = 20, type = "univariate", family = "gaussian")

  # stops at detection by default
  expect_equal(length(res$stat), res$detection_time)
})

test_that("sequential Gaussian: threshold=Inf runs to end of series", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_sequential_test(Y, threshold = Inf, type = "univariate", family = "gaussian")

  expect_equal(length(res$stat), length(Y))   # 1000
  expect_equal(res$n,            length(Y))
  expect_null(res$detection_time)
})

test_that("sequential Gaussian: statistic is zero at t = 1", {
  set.seed(1)
  res <- focus_sequential_test(rnorm(100), threshold = Inf,
                       type = "univariate", family = "gaussian")

  expect_equal(res$stat[1], 0)
})

test_that("sequential Gaussian: statistic is non-negative everywhere", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_sequential_test(Y, threshold = Inf, type = "univariate", family = "gaussian")

  expect_true(all(res$stat >= 0))
})

test_that("sequential Gaussian: no detection returned when signal is absent", {
  set.seed(7)
  Y   <- rnorm(200)   # pure noise, threshold very high
  res <- focus_sequential_test(Y, threshold = 1e6, type = "univariate", family = "gaussian")

  expect_null(res$detection_time)
  expect_null(res$detected_changepoint)
})


# ===========================================================================
# focus_sequential_test — known vs unknown pre-change parameter
# ===========================================================================

test_that("known theta0 yields larger max statistic than unknown (seed 45)", {
  set.seed(45)
  Y <- c(rnorm(1000, 0), rnorm(500, -1))

  res_k <- focus_sequential_test(Y, threshold = Inf, type = "univariate",
                         family = "gaussian", theta0 = 0)
  res_u <- focus_sequential_test(Y, threshold = Inf, type = "univariate",
                         family = "gaussian")

  expect_gt(max(res_k$stat), max(res_u$stat))

  # hard values from reference run
  expect_equal(round(max(res_k$stat), 4), 479.5316)
  expect_equal(round(max(res_u$stat), 4), 307.2719)
})


# ===========================================================================
# focus_sequential_test — one-sided detection
# ===========================================================================

test_that("right-sided detects increase; left-sided does not (seed 789)", {
  set.seed(789)
  Y <- c(rnorm(800, 0), rnorm(400, 1.5))

  res_right <- focus_sequential_test(Y, threshold = 30, type = "univariate_one_sided",
                             family = "gaussian", side = "right")
  res_left  <- focus_sequential_test(Y, threshold = 30, type = "univariate_one_sided",
                             family = "gaussian", side = "left")

  expect_equal(res_right$detection_time, 816)
  expect_null(res_left$detection_time)
  expect_equal(res_right$type, "univariate_one_sided")
})


# ===========================================================================
# focus_sequential_test — exponential family models
# ===========================================================================

test_that("Poisson family runs and peaks near true changepoint (seed 101)", {
  set.seed(101)
  Y_p <- c(rpois(500, 2), rpois(500, 6))
  res <- focus_sequential_test(Y_p, threshold = Inf, type = "univariate", family = "poisson")

  expect_equal(res$family,       "poisson")
  expect_equal(length(res$stat), 1000)
  expect_gt(which.max(res$stat), 900)   # peak in second half
})

test_that("Bernoulli family runs and returns correct length and gets the right stat (seed 123)", {
  set.seed(123)
  Y_b <- c(rbinom(500, 1, 0.2), rbinom(500, 1, 0.5))
  res <- focus_sequential_test(Y_b, threshold = Inf, type = "univariate", family = "bernoulli")

  expect_equal(res$family,       "bernoulli")
  expect_equal(length(res$stat), 1000)
  expect_true(all(res$stat >= 0))
  expect_equal(res$stat[1000],       63.855042146248138)
  expect_equal(res$stat[499],       1.521954710361967)
})


test_that("Gamma family stores shape in result (seed 124)", {
  set.seed(124)
  Y_g <- c(rgamma(500, shape = 2, scale = 2), rgamma(500, shape = 2, scale = 0.5))
  res <- focus_sequential_test(Y_g, threshold = Inf, type = "univariate",
                       family = "gamma", shape = 2, theta0 = 2)

  expect_equal(res$family,       "gamma")
  expect_equal(res$shape,        2)
  expect_equal(length(res$stat), 1000)
})

test_that("Gamma family errors when shape is not supplied", {
  set.seed(124)
  Y_g <- c(rgamma(200, shape = 2, scale = 2), rgamma(200, shape = 2, scale = 0.5))

  expect_error(
    focus_sequential_test(Y_g, threshold = 10, type = "univariate", family = "gamma")
  )
})


# ===========================================================================
# focus_sequential_test — NPFOCuS
# ===========================================================================

test_that("npfocus returns two-column stat matrix (seed 123)", {
  set.seed(123)
  Y_np   <- c(rnorm(500), rcauchy(100))
  quants <- qnorm(seq(0.01, 0.99, length.out = 5))

  res <- focus_sequential_test(Y_np, threshold = c(80, 25), type = "npfocus",
                       family = "npfocus", quantiles = quants)

  expect_equal(ncol(res$stat), 2)
  expect_equal(res$family,     "npfocus")
  expect_true(all(res$stat >= 0))
})


# ===========================================================================
# focus_sequential_test — ARP
# ===========================================================================

test_that("ARP detector detects mean shift in AR(2) series (seed 123)", {
  set.seed(123)
  ar_coefs <- c(0.7, -0.3)
  Y <- c(arima.sim(n = 300, model = list(ar = ar_coefs), sd = 1),
         2 + arima.sim(n = 300, model = list(ar = ar_coefs), sd = 1))

  res <- focus_sequential_test(Y, threshold = 20, family="arp", type = "arp",
                       rho = ar_coefs, mu0_arp = 0)

  expect_equal(res$type,            "arp")
  expect_equal(res$detection_time,   306)
  expect_equal(res$detected_changepoint, 297)
})

test_that("ARP statistics match a brute-force GLR for AR orders 1 to 3, known and unknown pre-change mean", {
  # Exact GLR for a change in mean of an AR(p) process with unknown pre-change
  # mean, computed over all changepoint locations. Aligned with focus_offline().
  arp_brute_force <- function(x, rho) {
    p <- length(rho)
    n <- length(x) - p
    y <- vapply(seq_len(n), function(i) x[i + p] - sum(rho * x[i + p - seq_len(p)]), numeric(1))
    v <- c(1, 1 - cumsum(rho), rep(1 - sum(rho), n - p - 1))
    vmax <- 1 - sum(rho)
    stat <- numeric(n - 1)
    coeffs <- matrix(0, nrow = n - 1, ncol = 6)
    for (i in seq_len(n)) {
      if (i < n) {
        for (k in i:(n - 1)) coeffs[k, ] <- coeffs[k, ] + c(vmax^2, 0, 0, -2 * vmax * y[i], 0, y[i]^2)
      }
      for (j in seq_len(p)) {
        if (i > j) {
          coeffs[i - j, ] <- coeffs[i - j, ] + c((vmax - v[j])^2, v[j]^2, 2 * v[j] * (vmax - v[j]),
                                                 -2 * (vmax - v[j]) * y[i], -2 * v[j] * y[i], y[i]^2)
        }
      }
      if (i > p + 1) {
        for (k in 1:(i - p - 1)) coeffs[k, ] <- coeffs[k, ] + c(0, vmax^2, 0, 0, -2 * vmax * y[i], y[i]^2)
      }
      if (i > 1) {
        lr0 <- sum(y[1:i]^2) - sum(y[1:i])^2 / i
        cf <- coeffs[1:(i - 1), , drop = FALSE]
        A <- cf[, 1]; B <- cf[, 2]; C <- cf[, 3]; D <- cf[, 4]; E <- cf[, 5]; FF <- cf[, 6]
        mu0 <- (C * E - 2 * B * D) / (4 * A * B - C^2)
        mu1 <- (C * D - 2 * A * E) / (4 * A * B - C^2)
        lr1 <- A * mu0^2 + B * mu1^2 + C * mu0 * mu1 + D * mu0 + E * mu1 + FF
        stat[i - 1] <- max(lr0 - lr1)
      }
    }
    c(rep(-1, p), stat)
  }

  # Exact GLR when the pre-change mean mu0 is known: after centring and
  # whitening, a change at tau gives weights 0 before the change, v[j] for the
  # j-th of the first p observations after it and 1 - sum(rho) afterwards.
  arp_brute_force_known <- function(x, rho, mu0) {
    p <- length(rho)
    n <- length(x) - p
    y <- vapply(seq_len(n), function(i) (x[i + p] - mu0) - sum(rho * (x[i + p - seq_len(p)] - mu0)), numeric(1))
    v <- c(1, 1 - cumsum(rho))[seq_len(p)]
    vmax <- 1 - sum(rho)
    stat <- numeric(n - 1)
    for (i in 2:n) {
      best <- -Inf
      for (tau in 1:(i - 1)) {
        j <- seq_len(i - tau)
        w <- ifelse(j <= p, v[pmin(j, p)], vmax)
        best <- max(best, sum(w * y[(tau + 1):i])^2 / sum(w^2))
      }
      stat[i - 1] <- best
    }
    c(rep(-1, p), stat)
  }

  for (rho in list(0.7, -0.5, c(0.8, -0.2), c(0.95, -0.1, 0.1))) {
    for (seed in 1:3) {
      set.seed(seed)
      Y <- rep(c(10, 35), each = 40) + arima.sim(list(ar = rho), n = 80)

      off <- as.vector(focus_offline(Y, threshold = Inf, type = "arp", rho = rho)$stat)
      det <- detector_create(type = "arp", rho = rho)
      on <- vapply(Y, function(y) {
        detector_update(det, y)
        get_statistics(det, family = "arp")$stat
      }, numeric(1))

      expect_equal(off, arp_brute_force(Y, rho), tolerance = 1e-6)
      expect_equal(on[-1], off)

      # known pre-change mean
      off_known <- as.vector(focus_offline(Y, threshold = Inf, type = "arp", rho = rho, mu0_arp = 10)$stat)
      det_known <- detector_create(type = "arp", rho = rho, mu0_arp = 10)
      on_known <- vapply(Y, function(y) {
        detector_update(det_known, y)
        get_statistics(det_known, family = "arp")$stat
      }, numeric(1))

      expect_equal(off_known, arp_brute_force_known(Y, rho, mu0 = 10), tolerance = 1e-6)
      expect_equal(on_known[-1], off_known)
    }
  }
})


# ===========================================================================
# Online mode vs offline
# ===========================================================================


test_that("sequential and offline produce identical stat traces", {
  set.seed(42)
  Y <- rnorm(200)

  res_seq <- focus_sequential_test(Y, threshold = Inf, type = "univariate", family = "gaussian")
  res_off <- focus_offline(Y, threshold = Inf, type = "univariate", family = "gaussian")

  expect_equal(as.vector(res_seq$stat), as.vector(res_off$stat))
})

test_that("online mode supports |> pipe chaining", {
  det <- detector_create(type = "univariate")
  set.seed(1)
  Y <- rnorm(5)

  det <- det |> detector_update(Y[1]) |> detector_update(Y[2]) |> detector_update(Y[3])
  expect_equal(detector_info_n(det), 3)
})


# ===========================================================================
# Detector inspection functions
# ===========================================================================

test_that("detector_info_n returns observation count", {
  det <- detector_create(type = "univariate")
  Y   <- make_gaussian(seed = 123)
  for (i in 1:10) detector_update(det, Y[i])

  expect_equal(detector_info_n(det), 10)
})

test_that("detector_info_sn returns cumulative sum (seed 123, first 10 obs)", {
  Y   <- make_gaussian(seed = 123)
  det <- detector_create(type = "univariate")
  for (i in 1:10) detector_update(det, Y[i])

  expect_equal(round(detector_info_sn(det), 6), 0.746256)
})

test_that("detector_cands_len returns candidate count (seed 123, first 10 obs)", {
  Y   <- make_gaussian(seed = 123)
  det <- detector_create(type = "univariate")
  for (i in 1:10) detector_update(det, Y[i])

  expect_equal(detector_cands_len(det), 7)
})

test_that("detector_candidates returns a data frame with correct row count", {
  Y   <- make_gaussian(seed = 123)
  det <- detector_create(type = "univariate")
  for (i in 1:10) detector_update(det, Y[i])

  cands <- detector_candidates(det)

  expect_s3_class(cands, "data.frame")
  expect_equal(nrow(cands), 7)
})

test_that("detector_cands_len grows with more observations", {
  det <- detector_create(type = "univariate")
  set.seed(1)
  Y <- rnorm(50)

  for (i in 1:10)  detector_update(det, Y[i])
  n10 <- detector_cands_len(det)

  for (i in 11:50) detector_update(det, Y[i])
  n50 <- detector_cands_len(det)

  # candidate set can grow (or at least not shrink to nothing)
  expect_gte(n10, 1)
  expect_gte(n50, 1)
})
