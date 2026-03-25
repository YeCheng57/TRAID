#' Extract genes showing ASE
#'
#' @param x An `ASEGeneResult` object.
#' @param sampleID Optional sample ID.
#'
#' @return A data.frame.
#' @export
getASEGenes <- function(x, sampleID = NULL) {
  if (!inherits(x, "ASEGeneResult")) {
    stop("`x` must be an ASEGeneResult object.")
  }

  df <- x$result

  if ("has_ASE" %in% colnames(df)) {
    df <- df[!is.na(df$has_ASE) & df$has_ASE, , drop = FALSE]
  }

  if (!is.null(sampleID)) {
    df <- df[df$sample_id == sampleID, , drop = FALSE]
  }

  df
}
