#' Methods for Detector, Statistics and Offline Result Objects
#'
#' Print, summary and plot methods for the S3 classes returned by the main
#' functions of the package: \code{"focus_detector"} objects, created by
#' \code{\link{detector_create}()}; \code{"focus_statistics"} objects, returned
#' by \code{\link{get_statistics}()}; and \code{"focus_offline"} objects,
#' returned by \code{\link{focus_offline}()}.
#'
#' @param x,object An object of class \code{"focus_detector"},
#'   \code{"focus_statistics"}, \code{"focus_offline"} or
#'   \code{"summary.focus_offline"}, as appropriate.
#' @param digits Number of significant digits used to print the statistics.
#' @param type,lty,col,xlab,ylab,main Graphical parameters passed to
#'   \code{\link[graphics]{matplot}()}. By default, one solid line is drawn
#'   per test statistic.
#' @param ... Further arguments passed to or from other methods (for
#'   \code{plot}, to \code{\link[graphics]{matplot}()}).
#'
#' @details
#' The \code{print} methods give a compact description of the object: for a
#' detector, its type, the number of observations processed and the number of
#' candidate changepoints currently stored; for the statistics, the current
#' time, changepoint estimate and test statistic(s); for an offline result, the
#' detector type and family, the threshold(s) and the detection, if any.
#'
#' The \code{summary} method for \code{"focus_offline"} objects additionally
#' reports, for each test statistic, its maximum over time, the time at which
#' the maximum is attained, the changepoint estimate at that time and the
#' threshold.
#'
#' The \code{plot} method for \code{"focus_offline"} objects draws the trace of
#' the test statistic(s) over time, with the finite threshold(s) as dashed
#' horizontal lines and, if a detection occurred, the detection time and the
#' estimated changepoint as dotted vertical lines.
#'
#' @return The \code{print} and \code{plot} methods return \code{x} invisibly.
#'   The \code{summary} method returns an object of class
#'   \code{"summary.focus_offline"}: a list with the elements \code{type},
#'   \code{family}, \code{shape}, \code{n}, \code{detection_time} and
#'   \code{detected_changepoint} of the offline result, the number of final
#'   candidates \code{n_candidates}, and \code{statistics}, a data frame with
#'   one row per test statistic.
#'
#' @examples
#' set.seed(123)
#' Y <- c(rnorm(100, mean = 0), rnorm(100, mean = 2))
#'
#' # Online interface
#' det <- detector_create(type = "univariate")
#' for (y in Y[1:120]) detector_update(det, y)
#' det
#' get_statistics(det, family = "gaussian")
#'
#' # Offline interface
#' res <- focus_offline(Y, threshold = 20, type = "univariate",
#'                      family = "gaussian")
#' res
#' summary(res)
#' plot(res)
#'
#' @name focus-methods
NULL

# Prints one "label: value" line, with the values aligned.
.focus_field <- function(label, value) {
  cat("  ", formatC(paste0(label, ":"), width = 15L, flag = "-"), value, "\n",
      sep = "")
}

# Formats each number separately, with its own significant digits.
.focus_format <- function(x, digits = 4L) {
  vapply(x, format, character(1L), digits = digits)
}

.focus_stat_names <- function(n_stats, family) {
  if (identical(family, "npfocus") && n_stats == 2L) {
    c("sum", "max")
  } else if (n_stats == 1L) {
    "stat"
  } else {
    paste0("stat", seq_len(n_stats))
  }
}

.focus_family_label <- function(family, shape) {
  if (identical(family, "gamma") && !is.null(shape)) {
    paste0("gamma (shape = ", .focus_format(shape), ")")
  } else {
    family
  }
}

.focus_detection_label <- function(detection_time, changepoint) {
  if (is.null(detection_time)) return("none")
  out <- paste0("at time ", .focus_format(detection_time))
  if (!is.null(changepoint)) {
    out <- paste0(out, " (changepoint estimate: ", .focus_format(changepoint), ")")
  }
  out
}

#' @rdname focus-methods
#' @export
print.focus_detector <- function(x, ...) {
  cat("focus detector\n")
  n <- tryCatch(detector_info_n(x), error = function(e) NULL)
  if (is.null(n)) {
    cat("  invalid: detectors cannot be restored from a saved session\n")
    return(invisible(x))
  }
  type <- attr(x, "type")
  side <- attr(x, "side")
  .focus_field("type", if (is.null(side)) type else paste0(type, " (side = \"", side, "\")"))
  .focus_field("observations", n)
  if (identical(type, "multivariate") && n > 0L) {
    .focus_field("dimensions", length(detector_info_sn(x)))
  }
  if (identical(type, "npfocus")) {
    .focus_field("quantiles", length(detector_info_sn(x)))
  }
  if (identical(type, "arp")) {
    .focus_field("AR order", attr(x, "ar_order"))
  } else {
    .focus_field("candidates", detector_cands_len(x))
  }
  invisible(x)
}

#' @rdname focus-methods
#' @export
print.focus_statistics <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  family <- attr(x, "family")
  cat("focus statistics", if (!is.null(family)) paste0(" (family: ", family, ")"),
      "\n", sep = "")
  .focus_field("stopping time", .focus_format(x$stopping_time, digits))
  .focus_field("changepoint", if (is.null(x$changepoint)) "not available"
               else .focus_format(x$changepoint, digits))
  if (is.null(x$stat)) {
    .focus_field("statistic", "not available")
  } else if (length(x$stat) == 1L) {
    .focus_field("statistic", .focus_format(x$stat, digits))
  } else {
    .focus_field("statistics", paste0(.focus_stat_names(length(x$stat), family), " = ",
                                      .focus_format(x$stat, digits), collapse = ", "))
  }
  invisible(x)
}

#' @rdname focus-methods
#' @export
print.focus_offline <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("focus offline detection\n")
  .focus_field("detector type", x$type)
  .focus_field("family", .focus_family_label(x$family, x$shape))
  .focus_field("observations", x$n)
  .focus_field("threshold", paste(.focus_format(x$threshold, digits), collapse = ", "))
  .focus_field("detection", .focus_detection_label(x$detection_time, x$detected_changepoint))
  invisible(x)
}

#' @rdname focus-methods
#' @export
summary.focus_offline <- function(object, ...) {
  stat <- as.matrix(object$stat)
  n_stats <- ncol(stat)
  if (nrow(stat) > 0L) {
    time_of_max <- max.col(t(stat), ties.method = "first")
    stat_max <- stat[cbind(time_of_max, seq_len(n_stats))]
    changepoint_at_max <- object$changepoint[time_of_max]
  } else {
    time_of_max <- changepoint_at_max <- rep(NA_integer_, n_stats)
    stat_max <- rep(NA_real_, n_stats)
  }
  statistics <- data.frame(
    max = stat_max,
    `time of max` = time_of_max,
    `changepoint at max` = changepoint_at_max,
    threshold = rep_len(object$threshold, n_stats),
    row.names = .focus_stat_names(n_stats, object$family),
    check.names = FALSE
  )
  structure(
    list(type = object$type, family = object$family, shape = object$shape,
         n = object$n, detection_time = object$detection_time,
         detected_changepoint = object$detected_changepoint,
         n_candidates = length(object$candidates$tau),
         statistics = statistics),
    class = "summary.focus_offline"
  )
}

#' @rdname focus-methods
#' @export
print.summary.focus_offline <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("focus offline detection: summary\n")
  .focus_field("detector type", x$type)
  .focus_field("family", .focus_family_label(x$family, x$shape))
  .focus_field("observations", x$n)
  .focus_field("detection", .focus_detection_label(x$detection_time, x$detected_changepoint))
  .focus_field("candidates", x$n_candidates)
  cat("\nStatistics:\n")
  print(x$statistics, digits = digits)
  invisible(x)
}

#' @rdname focus-methods
#' @export
plot.focus_offline <- function(x, type = "l", lty = 1, col = NULL, xlab = "Time",
                               ylab = "Statistic", main = NULL, ...) {
  stat <- as.matrix(x$stat)
  n_stats <- ncol(stat)
  if (is.null(col)) col <- seq_len(n_stats)
  graphics::matplot(seq_len(nrow(stat)), stat, type = type, lty = lty, col = col,
                    xlab = xlab, ylab = ylab, main = main, ...)
  threshold <- rep_len(x$threshold, n_stats)
  finite <- is.finite(threshold)
  if (length(x$threshold) == 1L) {
    if (finite[1L]) graphics::abline(h = x$threshold, lty = 2)
  } else if (any(finite)) {
    graphics::abline(h = threshold[finite], lty = 2, col = rep_len(col, n_stats)[finite])
  }
  if (!is.null(x$detected_changepoint)) {
    graphics::abline(v = x$detected_changepoint, lty = 3, col = "grey50")
  }
  if (n_stats > 1L) {
    graphics::legend("topleft", legend = .focus_stat_names(n_stats, x$family),
                     col = col, lty = lty, bty = "n")
  }
  invisible(x)
}

