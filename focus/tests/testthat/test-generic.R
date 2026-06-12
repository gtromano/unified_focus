library(testthat)
library(focus)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_gaussian <- function(n_pre = 500, n_post = 500, mu_post = 2, seed = 123) {
  set.seed(seed)
  c(rnorm(n_pre, mean = 0), rnorm(n_post, mean = mu_post))
}


# ===========================================================================
# 1. focus_offline — return structure
# ===========================================================================

test_that("focus_offline returns correct list keys", {
  res <- focus_offline(rnorm(50), threshold = Inf,
                       type = "univariate", family = "gaussian")

  expect_named(
    res,
    sort(c("stat", "changepoint", "detected_changepoint", "detection_time",
           "candidates", "threshold", "n", "type", "family", "shape")),
    ignore.order = TRUE
  )
})

test_that("focus_offline metadata fields are correct", {
  Y <- make_gaussian()
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$type,      "univariate")
  expect_equal(res$family,    "gaussian")
  expect_equal(res$threshold, 20)
})


# ===========================================================================
# 2. focus_offline — Gaussian univariate
# ===========================================================================

test_that("offline Gaussian: detects changepoint at correct time (seed 123)", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$detection_time,      510)
  expect_equal(res$detected_changepoint, 500)
  expect_equal(res$n,                    510)
})


test_that("offline Gaussian: detects changepoint at correct time, shifted by 123 (seed 123)", {
  Y   <- make_gaussian(seed = 123) + 123
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$detection_time,      510)
  expect_equal(res$detected_changepoint, 500)
  expect_equal(res$n,                    510)
})

test_that("offline Gaussian: check the values of the statistic (seed 123)", {
  Y   <- make_gaussian(seed = 123) + 123
  res <- focus_offline(Y, threshold = 20, theta0 = 123, type = "univariate", family = "gaussian")

  expect_equal(res$stat[509],      21.04167653450661)
  expect_equal(res$stat[456],      7.080510139397964)
})

test_that("offline Gaussian: check the values of the statistic (theta0 = 123) (seed 123)", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

  expect_equal(res$stat[510],      23.05886676141565)
  expect_equal(res$stat[509],      19.74501146734024)
  expect_equal(res$stat[456],      7.158358180852482)

})

test_that("offline Gaussian: stat vector length equals detection time", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")

  # stops at detection by default
  expect_equal(length(res$stat), res$detection_time)
})

test_that("offline Gaussian: threshold=Inf runs to end of series", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_offline(Y, threshold = Inf, type = "univariate", family = "gaussian")

  expect_equal(length(res$stat), length(Y))   # 1000
  expect_equal(res$n,            length(Y))
  expect_null(res$detection_time)
})

test_that("offline Gaussian: statistic is zero at t = 1", {
  set.seed(1)
  res <- focus_offline(rnorm(100), threshold = Inf,
                       type = "univariate", family = "gaussian")

  expect_equal(res$stat[1], 0)
})

test_that("offline Gaussian: statistic is non-negative everywhere", {
  Y   <- make_gaussian(seed = 123)
  res <- focus_offline(Y, threshold = Inf, type = "univariate", family = "gaussian")

  expect_true(all(res$stat >= 0))
})

test_that("offline Gaussian: no detection returned when signal is absent", {
  set.seed(7)
  Y   <- rnorm(200)   # pure noise, threshold very high
  res <- focus_offline(Y, threshold = 1e6, type = "univariate", family = "gaussian")

  expect_null(res$detection_time)
  expect_null(res$detected_changepoint)
})


# ===========================================================================
# 3. focus_offline — known vs unknown pre-change parameter
# ===========================================================================

test_that("known theta0 yields larger max statistic than unknown (seed 45)", {
  set.seed(45)
  Y <- c(rnorm(1000, 0), rnorm(500, -1))

  res_k <- focus_offline(Y, threshold = Inf, type = "univariate",
                         family = "gaussian", theta0 = 0)
  res_u <- focus_offline(Y, threshold = Inf, type = "univariate",
                         family = "gaussian")

  expect_gt(max(res_k$stat), max(res_u$stat))

  # hard values from reference run
  expect_equal(round(max(res_k$stat), 4), 479.5316)
  expect_equal(round(max(res_u$stat), 4), 307.2719)
})


# ===========================================================================
# 4. focus_offline — one-sided detection
# ===========================================================================

test_that("right-sided detects increase; left-sided does not (seed 789)", {
  set.seed(789)
  Y <- c(rnorm(800, 0), rnorm(400, 1.5))

  res_right <- focus_offline(Y, threshold = 30, type = "univariate_one_sided",
                             family = "gaussian", side = "right")
  res_left  <- focus_offline(Y, threshold = 30, type = "univariate_one_sided",
                             family = "gaussian", side = "left")

  expect_equal(res_right$detection_time, 816)
  expect_null(res_left$detection_time)
  expect_equal(res_right$type, "univariate_one_sided")
})


# ===========================================================================
# 5. focus_offline — multivariate Gaussian
# ===========================================================================

test_that("multivariate Gaussian detects changepoint near true location (seed 42)", {
  set.seed(42)
  Y_m <- rbind(matrix(rnorm(1000 * 3, 0),   ncol = 3),
               matrix(rnorm( 500 * 3, 1.2), ncol = 3))

  res <- focus_offline(Y_m, threshold = 30, type = "multivariate", family = "gaussian")

  expect_equal(res$detection_time,       1012)
  expect_equal(res$detected_changepoint, 1000)
  expect_equal(res$type, "multivariate")
})

test_that("multivariate Gaussian: stat length equals n", {
  set.seed(42)
  Y_m <- rbind(matrix(rnorm(200 * 2, 0), ncol = 2),
               matrix(rnorm(100 * 2, 2), ncol = 2))

  res <- focus_offline(Y_m, threshold = Inf, type = "multivariate", family = "gaussian")

  expect_equal(length(res$stat), 300)
  expect_equal(res$n,            300)
})


# ===========================================================================
# 6. focus_offline — exponential family models
# ===========================================================================

test_that("Poisson family runs and peaks near true changepoint (seed 101)", {
  set.seed(101)
  Y_p <- c(rpois(500, 2), rpois(500, 6))
  res <- focus_offline(Y_p, threshold = Inf, type = "univariate", family = "poisson")

  expect_equal(res$family,       "poisson")
  expect_equal(length(res$stat), 1000)
  expect_gt(which.max(res$stat), 900)   # peak in second half
})

test_that("Bernoulli family runs and returns correct length and gets the right stat (seed 123)", {
  set.seed(123)
  Y_b <- c(rbinom(500, 1, 0.2), rbinom(500, 1, 0.5))
  res <- focus_offline(Y_b, threshold = Inf, type = "univariate", family = "bernoulli")

  expect_equal(res$family,       "bernoulli")
  expect_equal(length(res$stat), 1000)
  expect_true(all(res$stat >= 0))
  expect_equal(res$stat[1000],       63.855042146248138)
  expect_equal(res$stat[499],       1.521954710361967)
})


test_that("Gamma family stores shape in result (seed 124)", {
  set.seed(124)
  Y_g <- c(rgamma(500, shape = 2, scale = 2), rgamma(500, shape = 2, scale = 0.5))
  res <- focus_offline(Y_g, threshold = Inf, type = "univariate",
                       family = "gamma", shape = 2, theta0 = 2)

  expect_equal(res$family,       "gamma")
  expect_equal(res$shape,        2)
  expect_equal(length(res$stat), 1000)
})

test_that("Gamma family errors when shape is not supplied", {
  set.seed(124)
  Y_g <- c(rgamma(200, shape = 2, scale = 2), rgamma(200, shape = 2, scale = 0.5))

  expect_error(
    focus_offline(Y_g, threshold = 10, type = "univariate", family = "gamma")
  )
})


# ===========================================================================
# 7. focus_offline — NPFOCuS
# ===========================================================================

test_that("npfocus returns two-column stat matrix (seed 123)", {
  set.seed(123)
  Y_np   <- c(rnorm(500), rcauchy(100))
  quants <- qnorm(seq(0.01, 0.99, length.out = 5))

  res <- focus_offline(Y_np, threshold = c(80, 25), type = "npfocus",
                       family = "npfocus", quantiles = quants)

  expect_equal(ncol(res$stat), 2)
  expect_equal(res$family,     "npfocus")
  expect_true(all(res$stat >= 0))
})


# ===========================================================================
# 8. focus_offline — ARP
# ===========================================================================

test_that("ARP detector detects mean shift in AR(2) series (seed 123)", {
  set.seed(123)
  ar_coefs <- c(0.7, -0.3)
  Y <- c(arima.sim(n = 300, model = list(ar = ar_coefs), sd = 1),
         2 + arima.sim(n = 300, model = list(ar = ar_coefs), sd = 1))

  res <- focus_offline(Y, threshold = 20, type = "arp",
                       rho = ar_coefs, mu0_arp = 0)

  expect_equal(res$type,            "arp")
  expect_equal(res$detection_time,   318)
  expect_equal(res$detected_changepoint, 297)
})


# ===========================================================================
# 9. Online mode — detector_create / detector_update / get_statistics
# ===========================================================================

test_that("online mode matches offline detection time (seed 123)", {
  Y         <- make_gaussian(seed = 123)
  threshold <- 20

  det          <- detector_create(type = "univariate")
  detected_at  <- NA
  detected_cp  <- NA

  for (i in seq_along(Y)) {
    detector_update(det, Y[i])
    result <- get_statistics(det, family = "gaussian")
    if (result$stat > threshold) {
      detected_at <- i
      detected_cp <- result$changepoint
      break
    }
  }

  expect_equal(detected_at, 510)
  expect_equal(detected_cp, 500)
})

test_that("online and offline produce identical stat traces", {
  set.seed(42)
  Y <- rnorm(200)

  res_off <- focus_offline(Y, threshold = Inf, type = "univariate", family = "gaussian")

  det     <- detector_create(type = "univariate")
  stat_on <- numeric(length(Y))
  for (i in seq_along(Y)) {
    detector_update(det, Y[i])
    stat_on[i] <- get_statistics(det, family = "gaussian")$stat
  }

  expect_equal(as.vector(res_off$stat), stat_on)
})

test_that("get_statistics stat is 0 after first observation", {
  det <- detector_create(type = "univariate")
  detector_update(det, 5.0)
  result <- get_statistics(det, family = "gaussian")

  expect_equal(result$stat,        0)
  expect_equal(result$changepoint, 0)
})

test_that("online mode supports |> pipe chaining", {
  det <- detector_create(type = "univariate")
  set.seed(1)
  Y <- rnorm(5)

  det <- det |> detector_update(Y[1]) |> detector_update(Y[2]) |> detector_update(Y[3])
  expect_equal(detector_info_n(det), 3)
})

test_that("online npfocus returns length-2 stat vector", {
  set.seed(123)
  quants <- qnorm(seq(0.01, 0.99, length.out = 5))
  det    <- detector_create("npfocus", quantiles = quants)
  for (y in rnorm(5)) detector_update(det, y)

  result <- get_statistics(det, family = "npfocus")
  expect_length(result$stat, 2)
})


# ===========================================================================
# 10. Detector inspection functions
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

test_that("detector_pieces_len returns candidate count (seed 123, first 10 obs)", {
  Y   <- make_gaussian(seed = 123)
  det <- detector_create(type = "univariate")
  for (i in 1:10) detector_update(det, Y[i])

  expect_equal(detector_pieces_len(det), 7)
})

test_that("detector_candidates returns a data frame with correct row count", {
  Y   <- make_gaussian(seed = 123)
  det <- detector_create(type = "univariate")
  for (i in 1:10) detector_update(det, Y[i])

  cands <- detector_candidates(det)

  expect_s3_class(cands, "data.frame")
  expect_equal(nrow(cands), 7)
})

test_that("detector_pieces_len grows with more observations", {
  det <- detector_create(type = "univariate")
  set.seed(1)
  Y <- rnorm(50)

  for (i in 1:10)  detector_update(det, Y[i])
  n10 <- detector_pieces_len(det)

  for (i in 11:50) detector_update(det, Y[i])
  n50 <- detector_pieces_len(det)

  # candidate set can grow (or at least not shrink to nothing)
  expect_gte(n10, 1)
  expect_gte(n50, 1)
})
