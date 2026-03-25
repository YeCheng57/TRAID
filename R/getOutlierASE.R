#' Extract ASE outliers
#'
#' @param x An `ASEResult` object.
#' @param sampleID Optional sample ID.
#'
#' @return A data.frame.
#' @export
getOutlierASE <- function(x, sampleID = NULL) {
  if (!inherits(x, "ASEResult")) {
    stop("`x` must be an ASEResult object.")
  }

  df <- x$result
  if ("is_outlier" %in% colnames(df)) {
    df <- df[!is.na(df$is_outlier) & df$is_outlier, , drop = FALSE]
  }

  if (!is.null(sampleID)) {
    df <- df[df$sample_id == sampleID, , drop = FALSE]
  }

  df
}
