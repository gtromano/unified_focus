#' Methods for Detector, Statistics, Offline Result and Projection Objects
#'
#' Methods for the S3 classes returned by the main functions of the package:
#' \code{"focus_detector"} objects, created by \code{\link{detector_create}()};
#' \code{"focus_statistics"} objects, returned by \code{\link{get_statistics}()};
#' \code{"focus_offline"} objects, returned by \code{\link{focus_offline}()};
#' and \code{"focus_projections"} objects, returned by
#' \code{\link{generate_projection_indexes}()}.
#'
#' @param x,object An object of class \code{"focus_detector"},
#'   \code{"focus_statistics"}, \code{"focus_offline"},
#'   \code{"focus_projections"}, \code{"summary.focus_detector"} or
#'   \code{"summary.focus_offline"}, as appropriate.
#' @param digits Number of significant digits used to print the statistics.
#' @param family,theta0,shape Optional arguments passed to
#'   \code{\link{get_statistics}()}: if \code{family} is given, the summary of
#'   a detector includes its current statistics.
#' @param data Optionally, the data passed to \code{\link{focus_offline}()},
#'   drawn above the trace of the statistic(s).
#' @param type,lty,col,xlab,ylab,main Graphical parameters passed to
#'   \code{\link[graphics]{matplot}()}. By default, one solid line is drawn
#'   per test statistic.
#' @param i Indices of the projections to keep.
#' @param n Maximum number of projections to print.
#' @param ... Further arguments passed to or from other methods (for
#'   \code{plot}, to \code{\link[graphics]{matplot}()}).
#'
#' @details
#' The \code{print} methods give a compact description of the object: for a
#' detector, its type, the number of observations processed and the number of
#' candidate changepoints currently stored; for the statistics, the current
#' time, changepoint estimate and test statistic(s); for an offline result, the
#' detector type and family, the threshold(s) and the detection, if any; for
#' projection index sets, the (0-based) indices of each projection, one
#' projection per row.
#'
#' The \code{summary} method for \code{"focus_detector"} objects additionally
#' reports the cumulative sums (see \code{\link{detector_info_sn}()}) and the
#' locations of the candidate changepoints, at the time it is called, and, if a
#' \code{family} is given, the current statistics.
#'
#' The \code{summary} method for \code{"focus_offline"} objects additionally
#' reports, for each test statistic, its maximum over time, the time at which
#' the maximum is attained, the changepoint estimate at that time and the
#' threshold.
#'
#' The \code{plot} method for \code{"focus_offline"} objects draws the trace of
#' the test statistic(s) over time, with the finite threshold(s) as dashed
#' horizontal lines and, if a detection occurred, the estimated changepoint as
#' a dotted vertical line. If \code{data} are given, they are drawn above the
#' trace, with one panel per dimension. If \pkg{ggplot2} is installed, \code{autoplot()} returns the same plot
#' as a \code{"ggplot"} object.
#'
#' Projection index sets can be subset with \code{[} (e.g., with
#' \code{head()}), keeping their class, and \code{as.matrix} returns them as a
#' matrix with one row per projection.
#'
#' @return The \code{print} and \code{plot} methods return \code{x} invisibly.
#'   The \code{summary} method for detectors returns an object of class
#'   \code{"summary.focus_detector"}, a list with the reported elements. The
#'   \code{summary} method for offline results returns an object of class
#'   \code{"summary.focus_offline"}: a list with the elements \code{type},
#'   \code{family}, \code{shape}, \code{n}, \code{detection_time} and
#'   \code{detected_changepoint} of the offline result, the number of final
#'   candidates \code{n_candidates}, and \code{statistics}, a data frame with
#'   one row per test statistic. The \code{[} method returns a
#'   \code{"focus_projections"} object and \code{as.matrix} an integer matrix.
#'
#' @examples
#' set.seed(123)
#' Y <- c(rnorm(100, mean = 0), rnorm(100, mean = 2))
#'
#' # Online interface
#' det <- detector_create(type = "univariate")
#' for (y in Y[1:120]) detector_update(det, y)
#' det
#' summary(det, family = "gaussian")
#' get_statistics(det, family = "gaussian")
#'
#' # Offline interface
#' res <- focus_offline(Y, threshold = 20, type = "univariate",
#'                      family = "gaussian")
#' res
#' summary(res)
#' plot(res)
#' plot(res, data = Y)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::autoplot(res, data = Y)
#' }
#'
#' # Projection index sets
#' proj <- generate_projection_indexes(d = 6, p = 2)
#' proj
#' head(proj, 2)
#' as.matrix(proj)
#'
#' @name focus-methods
NULL

# Prints one "label: value" line, with the values aligned.
.focus_field <- function(label, value) {
  cat("  ", formatC(paste0(label, ":"), width = 15L, flag = "-"), value, "\n",
      sep = "")
}

# Panel labels for the data: the column names, if any, or "Data" (univariate)
# and "Dimension j" (multivariate).
.focus_data_labels <- function(data) {
  if (!is.null(colnames(data))) return(colnames(data))
  if (ncol(data) == 1L) "Data" else paste("Dimension", seq_len(ncol(data)))
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
  }
  .focus_field("candidates", detector_cands_len(x))
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
plot.focus_offline <- function(x, data = NULL, type = "l", lty = 1, col = NULL, xlab = "Time",
                               ylab = "Statistic", main = NULL, ...) {
  stat <- as.matrix(x$stat)
  n_stats <- ncol(stat)
  if (is.null(col)) col <- seq_len(n_stats)
  xlim <- NULL
  if (!is.null(data)) {
    data <- as.matrix(data)
    labels <- .focus_data_labels(data)
    xlim <- c(1, max(nrow(data), nrow(stat)))
    # One panel per dimension of the data, above the trace of the statistic(s)
    # cex is set after mfrow, which would otherwise shrink the text with 3+ panels
    oldpar <- graphics::par(mfrow = c(ncol(data) + 1L, 1L),
                            mar = c(if (ncol(data) > 1L) 2 else 4, 4, 1, 1) + 0.1,
                            cex = if (ncol(data) > 1L) 0.8 else 1)
    on.exit(graphics::par(oldpar))
    for (j in seq_len(ncol(data))) {
      graphics::plot(seq_len(nrow(data)), data[, j], type = "l", col = "grey30", xlim = xlim,
                     xlab = "", ylab = labels[j], main = if (j == 1L) main)
      if (!is.null(x$detected_changepoint)) {
        graphics::abline(v = x$detected_changepoint, lty = 3, col = "grey50")
      }
    }
    graphics::par(mar = c(4, 4, 1, 1) + 0.1)
    main <- NULL
  }
  graphics::matplot(seq_len(nrow(stat)), stat, type = type, lty = lty, col = col, xlim = xlim,
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

# Formats a vector, printing at most `max` values.
.focus_format_values <- function(x, max = 8L) {
  x <- .focus_format(x)
  if (length(x) > max) x <- c(x[seq_len(max - 2L)], "...", x[length(x)])
  paste(x, collapse = ", ")
}

#' @rdname focus-methods
#' @export
summary.focus_detector <- function(object, family = NULL, theta0 = NULL, shape = NULL, ...) {
  n <- detector_info_n(object)
  type <- attr(object, "type")
  structure(
    list(type = type, side = attr(object, "side"), ar_order = attr(object, "ar_order"),
         n = n, sn = detector_info_sn(object),
         n_candidates = detector_cands_len(object),
         candidates = sort(unique(detector_candidates(object)$tau)),
         statistics = if (!is.null(family)) {
           get_statistics(object, family = family, theta0 = theta0, shape = shape)
         }),
    class = "summary.focus_detector"
  )
}

#' @rdname focus-methods
#' @export
print.summary.focus_detector <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("focus detector: summary\n")
  .focus_field("type", if (is.null(x$side)) x$type
               else paste0(x$type, " (side = \"", x$side, "\")"))
  .focus_field("observations", x$n)
  if (!is.null(x$ar_order)) .focus_field("AR order", x$ar_order)
  if (length(x$sn) > 0L) {
    .focus_field(if (length(x$sn) == 1L) "running sum" else "running sums",
                 .focus_format_values(signif(x$sn, digits)))
  }
  if (!is.null(x$n_candidates)) {
    .focus_field("candidates", x$n_candidates)
    if (length(x$candidates) > 0L) {
      .focus_field("locations", .focus_format_values(x$candidates))
    }
  }
  if (!is.null(x$statistics)) {
    cat("\n")
    print(x$statistics, digits = digits)
  }
  invisible(x)
}

#' @rdname focus-methods
#' @export
print.focus_projections <- function(x, n = 10L, ...) {
  cat("focus projection indexes (0-based)\n")
  .focus_field("dimensions", attr(x, "d"))
  .focus_field("projections", paste(length(x), "of size", attr(x, "p")))
  shown <- seq_len(min(length(x), n))
  labels <- formatC(paste0("[", shown, "]"), width = nchar(length(x)) + 2L)
  for (k in shown) cat("  ", labels[k], " ", paste(x[[k]], collapse = " "), "\n", sep = "")
  if (length(x) > n) cat("  ... and", length(x) - n, "more\n")
  invisible(x)
}

#' @rdname focus-methods
#' @export
`[.focus_projections` <- function(x, i) {
  structure(unclass(x)[i], d = attr(x, "d"), p = attr(x, "p"), class = "focus_projections")
}

#' @rdname focus-methods
#' @export
as.matrix.focus_projections <- function(x, ...) {
  out <- do.call(rbind, unclass(x))
  if (is.null(out)) out <- matrix(integer(0), nrow = 0L, ncol = attr(x, "p"))
  out
}

#' @rdname focus-methods
#' @exportS3Method ggplot2::autoplot
autoplot.focus_offline <- function(object, data = NULL, ...) {
  .data <- ggplot2::.data
  labels <- if (!is.null(data)) .focus_data_labels(as.matrix(data))
  panels <- c(labels, "Statistic")
  stat <- as.matrix(object$stat)
  stat_names <- .focus_stat_names(ncol(stat), object$family)
  trace <- data.frame(time = rep(seq_len(nrow(stat)), ncol(stat)), value = as.vector(stat),
                      statistic = factor(rep(stat_names, each = nrow(stat)), levels = stat_names),
                      panel = factor("Statistic", levels = panels))
  threshold <- rep_len(object$threshold, ncol(stat))
  thresholds <- data.frame(statistic = factor(stat_names, levels = stat_names),
                           threshold = threshold,
                           panel = factor("Statistic", levels = panels))[is.finite(threshold), ,
                                                                          drop = FALSE]
  multiple <- ncol(stat) > 1L
  p <- ggplot2::ggplot(trace, ggplot2::aes(x = .data$time, y = .data$value))
  if (!is.null(data)) {
    data <- as.matrix(data)
    # One panel per dimension of the data, above the trace of the statistic(s)
    observed <- data.frame(time = rep(seq_len(nrow(data)), ncol(data)), value = as.vector(data),
                           panel = factor(rep(labels, each = nrow(data)), levels = panels))
    p <- p + ggplot2::geom_line(data = observed, colour = "grey30")
  }
  if (multiple) {
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = .data$statistic))
  } else {
    p <- p + ggplot2::geom_line()
  }
  if (nrow(thresholds) > 0L) {
    p <- p + if (multiple) {
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$threshold, colour = .data$statistic),
                          data = thresholds, linetype = "dashed")
    } else {
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$threshold), data = thresholds,
                          linetype = "dashed")
    }
  }
  if (!is.null(object$detected_changepoint)) {
    p <- p + ggplot2::geom_vline(xintercept = object$detected_changepoint, linetype = "dotted",
                                 colour = "grey50")
  }
  if (!is.null(data)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), ncol = 1, scales = "free_y")
  }
  p + ggplot2::labs(x = "Time", y = if (is.null(data)) "Statistic" else NULL, colour = NULL)
}
