#' Summarize ASResult object
#'
#' @param object An `ASResult` object.
#' @param ... Reserved.
#'
#' @return A summary object of class `summary.ASResult`.
#' @export
summary.ASResult <- function(object, ...) {
  if (!inherits(object, "ASResult")) {
    stop("`object` must be an ASResult.")
  }

  start_df <- object$result_start
  end_df <- object$result_end

  out <- list(
    mode = object$metadata$mode,
    returnAll = object$metadata$returnAll,
    n_samples = object$metadata$n_samples,
    n_input_rows = object$metadata$n_rows_input,
    n_start_rows = if (is.null(start_df)) 0L else nrow(start_df),
    n_end_rows = if (is.null(end_df)) 0L else nrow(end_df),
    n_start_outliers = if (is.null(start_df)) 0L else sum(start_df$is_outlier, na.rm = TRUE),
    n_end_outliers = if (is.null(end_df)) 0L else sum(end_df$is_outlier, na.rm = TRUE)
  )

  class(out) <- "summary.ASResult"
  out
}

#' @export
print.summary.ASResult <- function(x, ...) {
  cat("ASResult summary\n")
  cat(" Mode           :", x$mode, "\n")
  cat(" returnAll      :", x$returnAll, "\n")
  cat(" Samples        :", x$n_samples, "\n")
  cat(" Input rows     :", x$n_input_rows, "\n")
  cat(" Start rows     :", x$n_start_rows, "\n")
  cat(" End rows       :", x$n_end_rows, "\n")
  cat(" Start outliers :", x$n_start_outliers, "\n")
  cat(" End outliers   :", x$n_end_outliers, "\n")
  invisible(x)
}
#' Plot ASResult object
#'
#' @param x An `ASResult` object.
#' @param type Type of plot: "zscore", "junction", "sample_summary".
#' @param site One of "start" or "end".
#' @param sampleID Required for zscore plot; optional for junction plot.
#' @param chrom,start,end Required for junction plot.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#' @param ... Additional arguments passed to underlying plotting functions.
#'
#' @return A ggplot object.
#' @export
plot.ASResult <- function(
    x,
    type = c("zscore", "junction", "sample_summary"),
    site = c("start", "end"),
    sampleID = NULL,
    chrom = NULL,
    start = NULL,
    end = NULL,
    groupBy = NULL,
    ...
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  type <- match.arg(type)
  site <- match.arg(site)

  if (type == "zscore") {
    if (is.null(sampleID)) {
      stop("`sampleID` is required for zscore plot.")
    }
    return(plotASZScore(x, sampleID = sampleID, site = site, ...))
  }

  if (type == "junction") {
    if (is.null(chrom) || is.null(start) || is.null(end)) {
      stop("`chrom`, `start`, and `end` are required for junction plot.")
    }
    return(plotASJunction(
      x,
      chrom = chrom,
      start = start,
      end = end,
      site = site,
      sampleID = sampleID,
      groupBy = groupBy,
      ...
    ))
  }

  if (type == "sample_summary") {
    return(plotASSampleSummary(x, site = site, groupBy = groupBy, ...))
  }

  stop("Unknown plot type.")
}
