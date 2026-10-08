library(testthat)
library(focus)

# ===========================================================================
# S3 classes and methods
# ===========================================================================

make_detector <- function(type = "univariate", Y = NULL, ...) {
  if (is.null(Y)) {
    set.seed(1)
    Y <- c(rnorm(100), rnorm(50, mean = 1))
  }
  det <- detector_create(type = type, ...)
  for (i in seq_len(NROW(Y))) detector_update(det, if (is.matrix(Y)) Y[i, ] else Y[i])
  det
}

# Runs a plot method on a null device.
plot_quietly <- function(...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  plot(...)
}

test_that("the inspection functions are S3 methods for focus_detector objects", {
  generics <- attr(methods(class = "focus_detector"), "info")$generic
  expect_setequal(generics, c("detector_candidates", "detector_cands_len", "detector_info_n",
                              "detector_info_sn", "print", "summary"))
  expect_error(detector_info_n(list()), "no applicable method")
  # detector_update() and get_statistics() are plain functions, checking their input
  expect_error(get_statistics(list(), family = "gaussian"), "must be a detector")
  expect_error(detector_update(1, 2), "must be a detector")
})

test_that("detector_update returns the detector, updated in place", {
  det <- detector_create(type = "univariate")
  out <- detector_update(det, 0.5)
  expect_s3_class(out, "focus_detector")
  expect_identical(out, det)
  detector_update(det, 0.1, lambda = 2)
  expect_equal(detector_info_n(det), 3)
})

test_that("get_statistics returns a focus_statistics object", {
  det <- make_detector()
  r <- get_statistics(det, family = "gaussian")
  expect_s3_class(r, "focus_statistics")
  expect_equal(attr(r, "family"), "gaussian")
  expect_named(r, c("stopping_time", "changepoint", "stat"))
  expect_output(print(r), "focus statistics \\(family: gaussian\\)")
})

test_that("summary of a detector reports its state", {
  det <- make_detector()
  s <- summary(det)
  expect_s3_class(s, "summary.focus_detector")
  expect_equal(s$n, 150)
  expect_equal(s$sn, detector_info_sn(det))
  expect_equal(s$n_candidates, detector_cands_len(det))
  expect_equal(s$candidates, sort(unique(detector_candidates(det)$tau)))
  expect_null(s$statistics)
  expect_output(print(s), "focus detector: summary")

  s2 <- summary(det, family = "gaussian")
  expect_equal(s2$statistics, get_statistics(det, family = "gaussian"))
  expect_output(print(s2), "focus statistics")

  det_arp <- make_detector("arp", rho = 0.5)
  s_arp <- summary(det_arp, family = "arp")
  expect_equal(s_arp$ar_order, 1L)
  expect_equal(s_arp$n_candidates, detector_cands_len(det_arp))
  expect_gt(s_arp$n_candidates, 0L)
  expect_equal(s_arp$candidates, sort(unique(detector_candidates(det_arp)$tau)))
  expect_output(print(s_arp), "AR order")
  expect_output(print(det_arp), "candidates")
})

test_that("offline results have print, summary and plot methods", {
  set.seed(123)
  Y <- c(rnorm(100), rnorm(100, mean = 2))
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")
  expect_s3_class(res, "focus_offline")
  expect_output(print(res), "focus offline detection")

  s <- summary(res)
  expect_s3_class(s, "summary.focus_offline")
  expect_equal(s$statistics$max, max(res$stat))
  expect_output(print(s), "Statistics:")

  expect_invisible(plot_quietly(res))
  expect_silent(plot_quietly(res, data = Y))

  res_np <- focus_offline(Y, threshold = c(10, 5), type = "npfocus", family = "npfocus",
                          quantiles = qnorm(c(0.25, 0.5, 0.75)))
  expect_silent(plot_quietly(res_np, data = Y))

  # Multivariate data: one panel per dimension
  set.seed(1)
  Y_mv <- rbind(matrix(rnorm(300), ncol = 3), matrix(rnorm(150, mean = 1), ncol = 3))
  res_mv <- focus_offline(Y_mv, threshold = 30, type = "multivariate")
  expect_silent(plot_quietly(res_mv, data = Y_mv))
})

test_that("generate_projection_indexes returns a focus_projections object", {
  proj <- generate_projection_indexes(6, 2)
  expect_s3_class(proj, "focus_projections")
  expect_length(proj, 6)
  expect_equal(attr(proj, "d"), 6L)
  expect_equal(attr(proj, "p"), 2L)
  expect_equal(proj[[6]], c(5L, 0L))
  expect_output(print(proj), "6 of size 2")

  first <- head(proj, 2)
  expect_s3_class(first, "focus_projections")
  expect_length(first, 2)
  expect_equal(as.matrix(proj), rbind(0:5, c(1:5, 0L)) |> t())
  expect_equal(dim(as.matrix(first)), c(2L, 2L))

  set.seed(4)
  Y <- matrix(rnorm(600), ncol = 6)
  det_cls <- make_detector("multivariate", Y = Y, dim_indexes = proj)
  det_lst <- make_detector("multivariate", Y = Y, dim_indexes = unclass(proj))
  expect_equal(get_statistics(det_cls, family = "gaussian"),
               get_statistics(det_lst, family = "gaussian"))
  det <- detector_create(type = "multivariate", dim_indexes = first)
  expect_s3_class(detector_update(det, rnorm(6)), "focus_detector")
})

test_that("autoplot returns a ggplot of an offline result", {
  skip_if_not_installed("ggplot2")
  set.seed(123)
  Y <- c(rnorm(100), rnorm(100, mean = 2))
  res <- focus_offline(Y, threshold = 20, type = "univariate", family = "gaussian")
  p <- ggplot2::autoplot(res)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  p_data <- ggplot2::autoplot(res, data = Y)
  expect_s3_class(p_data, "ggplot")
  expect_silent(ggplot2::ggplot_build(p_data))

  res_np <- focus_offline(Y, threshold = c(10, 5), type = "npfocus", family = "npfocus",
                          quantiles = qnorm(c(0.25, 0.5, 0.75)))
  expect_silent(ggplot2::ggplot_build(ggplot2::autoplot(res_np, data = Y)))

  # Multivariate data: one panel per dimension
  set.seed(1)
  Y_mv <- rbind(matrix(rnorm(300), ncol = 3), matrix(rnorm(150, mean = 1), ncol = 3))
  res_mv <- focus_offline(Y_mv, threshold = 30, type = "multivariate")
  built <- ggplot2::ggplot_build(ggplot2::autoplot(res_mv, data = Y_mv))
  expect_equal(nlevels(built$layout$layout$PANEL), 4L)
})
