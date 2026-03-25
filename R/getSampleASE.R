#' Extract ASE results for one sample
#'
#' @param x An `ASEResult` object.
#' @param sampleID Sample ID.
#'
#' @return A data.frame.
#' @export
getSampleASE <- function(x, sampleID) {
  if (!inherits(x, "ASEResult")) {
    stop("`x` must be an ASEResult object.")
  }
  x$result[x$result$sample_id == sampleID, , drop = FALSE]
}
