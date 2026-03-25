#' Extract gene-level ASE results for one gene
#'
#' @param x An `ASEGeneResult` object.
#' @param geneID Gene identifier or gene symbol.
#'
#' @return A data.frame.
#' @export
getGeneASEGene <- function(x, geneID) {
  if (!inherits(x, "ASEGeneResult")) {
    stop("`x` must be an ASEGeneResult object.")
  }

  df <- x$result
  if (!"gene_label" %in% colnames(df)) {
    stop("Column `gene_label` not found in ASEGeneResult.")
  }

  df[df$gene_label == geneID, , drop = FALSE]
}
