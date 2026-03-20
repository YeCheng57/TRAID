#' Extract AE results for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Character scalar specifying the target sample.
#'
#' @return A data.frame containing AE results for one sample.
#' @export
getSampleAE <- function(x, sampleID) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  if (!"sample_id" %in% colnames(x$result)) {
    stop("No `sample_id` column found in result.")
  }

  if (!sampleID %in% x$result$sample_id) {
    stop("`sampleID` not found in result.")
  }

  x$result[x$result$sample_id == sampleID, , drop = FALSE]
}
