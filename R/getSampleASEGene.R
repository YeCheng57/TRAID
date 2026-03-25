#' Extract gene-level ASE results for one sample
#'
#' @param x An `ASEGeneResult` object.
#' @param sampleID Sample ID.
#'
#' @return A data.frame.
#' @export
getSampleASEGene <- function(x, sampleID) {
  if (!inherits(x, "ASEGeneResult")) {
    stop("`x` must be an ASEGeneResult object.")
  }

  x$result[x$result$sample_id == sampleID, , drop = FALSE]
}
