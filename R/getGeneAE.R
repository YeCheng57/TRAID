#' Extract AE results for one gene
#'
#' @param x An `AEResult` object.
#' @param geneID Gene identifier or gene symbol.
#' @param use One of `"gene_id"` or `"gene_name"`.
#'
#' @return A data.frame.
#' @export
getGeneAE <- function(x, geneID, use = c("gene_id", "gene_name")) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  use <- match.arg(use)
  df <- x$result

  if (!use %in% colnames(df)) {
    stop("Column `", use, "` not found in AE result.")
  }

  df[df[[use]] == geneID, , drop = FALSE]
}
