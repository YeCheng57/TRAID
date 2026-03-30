#' Extract AS outliers
#'
#' @param x An `ASResult` object.
#' @param site One of "start", "end", "both".
#' @param sampleID Optional sample ID.
#'
#' @return A data.frame if one site is requested, or a list if `site = "both"`.
#' @export
getOutlierAS <- function(x, site = c("both","start", "end"), sampleID = NULL) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)

  extract_one <- function(res, sampleID = NULL) {
    if (is.null(res)) {
      return(NULL)
    }

    if ("is_outlier" %in% colnames(res)) {
      res <- res[!is.na(res$is_outlier) & res$is_outlier, , drop = FALSE]
    }

    if (!is.null(sampleID)) {
      res <- res[res$sample_id == sampleID, , drop = FALSE]
    }

    res
  }

  if (site == "start") {
    return(extract_one(x$result_start, sampleID))
  }

  if (site == "end") {
    return(extract_one(x$result_end, sampleID))
  }

  rbind(
    start = extract_one(x$result_start, sampleID),
    end = extract_one(x$result_end, sampleID)
  )
}
