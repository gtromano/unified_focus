#' Get the Number of Candidate Segments
#'
#' Returns the number of candidate changepoint segments currently tracked
#' by the detector.
#'
#' @param det_ptr A \code{"focus_detector"} object created by
#'   \code{\link{detector_create}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Integer. Number of candidate segments.
#'
#' @details
#' The FOCuS algorithm maintains a set of candidate segments that could
#' potentially contain changepoints. This number grows with time but is
#' controlled by the pruning parameters. \code{detector_cands_len()} is an S3
#' generic, with a method for \code{"focus_detector"} objects.
#'
#' @examples
#' set.seed(1)
#' det <- detector_create(type = "univariate")
#' for (y in rnorm(50)) detector_update(det, y)
#' detector_cands_len(det)
#'
#' @export
detector_cands_len <- function(det_ptr, ...) UseMethod("detector_cands_len")

#' @rdname detector_cands_len
#' @export
detector_cands_len.focus_detector <- function(det_ptr, ...) {
  .detector_cands_len(det_ptr)
}

#' Get the Number of Observations Processed
#'
#' Returns the total number of observations processed by the detector.
#'
#' @param det_ptr A \code{"focus_detector"} object created by
#'   \code{\link{detector_create}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Integer. Number of observations processed (current time index).
#'
#' @details
#' \code{detector_info_n()} is an S3 generic, with a method for
#' \code{"focus_detector"} objects.
#'
#' @examples
#' det <- detector_create(type = "univariate")
#' for (y in c(0.5, 1.2, -0.3)) detector_update(det, y)
#' detector_info_n(det)
#'
#' @export
detector_info_n <- function(det_ptr, ...) UseMethod("detector_info_n")

#' @rdname detector_info_n
#' @export
detector_info_n.focus_detector <- function(det_ptr, ...) {
  .detector_info_n(det_ptr)
}

#' Get the Cumulative Sum Statistic
#'
#' Returns the current cumulative sum statistic maintained by the detector.
#'
#' @param det_ptr A \code{"focus_detector"} object created by
#'   \code{\link{detector_create}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Numeric vector. Cumulative sum statistic. For univariate detectors,
#'   a scalar (length-1 vector). For multivariate detectors, a vector of
#'   length equal to the number of dimensions.
#'
#' @details
#' \code{detector_info_sn()} is an S3 generic, with a method for
#' \code{"focus_detector"} objects.
#'
#' @examples
#' det <- detector_create(type = "multivariate")
#' detector_update(det, c(0.5, 1.2))
#' detector_update(det, c(-0.3, 0.4))
#' detector_info_sn(det)
#'
#' @export
detector_info_sn <- function(det_ptr, ...) UseMethod("detector_info_sn")

#' @rdname detector_info_sn
#' @export
detector_info_sn.focus_detector <- function(det_ptr, ...) {
  .detector_info_sn(det_ptr)
}

#' Get the Candidate Segments
#'
#' Returns detailed information about all candidate changepoint segments
#' currently tracked by the detector.
#'
#' @param det_ptr A \code{"focus_detector"} object created by
#'   \code{\link{detector_create}()}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A data frame (tibble) with one row per candidate and columns:
#'   \item{tau}{Numeric vector. Candidate changepoint locations, on the same
#'     scale as the changepoint estimate returned by
#'     \code{\link{get_statistics}()}.}
#'   \item{st}{List of numeric vectors. Sufficient statistics for each
#'     candidate segment (e.g., cumulative sums of the data).}
#'   \item{side}{Character vector. Side indicator for each candidate
#'     (relevant for one-sided detectors).}
#'
#' @details
#' Each row represents a candidate segment from time \code{tau} to the current
#' time. The sufficient statistics in \code{st} are used to efficiently compute
#' test statistics without reprocessing the data. \code{detector_candidates()}
#' is an S3 generic, with a method for \code{"focus_detector"} objects.
#'
#' @examples
#' set.seed(1)
#' det <- detector_create(type = "univariate")
#' for (y in rnorm(50)) detector_update(det, y)
#' head(detector_candidates(det))
#'
#' @export
detector_candidates <- function(det_ptr, ...) UseMethod("detector_candidates")

#' @rdname detector_candidates
#' @export
detector_candidates.focus_detector <- function(det_ptr, ...) {
  .detector_candidates(det_ptr)
}
