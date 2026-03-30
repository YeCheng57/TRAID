#' Extract AS results for one sample
#'
#' @param x An `ASResult` object.
#' @param sampleID Character scalar specifying the target sample.
#' @param site One of "start", "end", "both".
#'
#' @return A data.frame if one site is requested, or a list if `site = "both"`.
#' @export
getSampleAS <- function(x, sampleID, site = c("start", "end", "both")) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)

  extract_one <- function(res, sampleID) {
    if (is.null(res)) {
      return(NULL)
    }
    if (!sampleID %in% res$sample_id) {
      return(res[0, , drop = FALSE])
    }
    res[res$sample_id == sampleID&res$is_outlier, , drop = FALSE]
  }

  if (site == "start") {
    return(extract_one(x$result_start, sampleID))
  }

  if (site == "end") {
    return(extract_one(x$result_end, sampleID))
  }

  list(
    start = extract_one(x$result_start, sampleID),
    end = extract_one(x$result_end, sampleID)
  )
}
