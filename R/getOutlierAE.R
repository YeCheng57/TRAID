#' Extract aberrant expression outliers
#'
#' @param x An `AEResult` object.
#' @param sampleID Optional character scalar specifying one sample.
#'
#' @return A data.frame containing outlier AE results.
#' @export
getOutlierAE <- function(x, sampleID = NULL) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  res <- x$result

  if (!"is_outlier" %in% colnames(res)) {
    stop("No `is_outlier` column found in result.")
  }

  res <- res[!is.na(res$is_outlier) & res$is_outlier, , drop = FALSE]

  if (!is.null(sampleID)) {
    if (!sampleID %in% res$sample_id) {
      stop("`sampleID` not found among outlier results.")
    }
    res <- res[res$sample_id == sampleID, , drop = FALSE]
  }

  res
}
