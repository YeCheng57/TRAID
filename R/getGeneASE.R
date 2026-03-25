#' Extract ASE site-level results for one gene
#'
#' @param x An `ASEResult` object.
#' @param geneID Gene identifier or gene symbol.
#' @param use One of `"gene_name"` or `"gene_id"`.
#'
#' @return A data.frame.
#' @export
getGeneASE <- function(x, geneID, use = c("gene_name", "gene_id")) {
  if (!inherits(x, "ASEResult")) {
    stop("`x` must be an ASEResult object.")
  }

  use <- match.arg(use)
  df <- x$result

  if (!use %in% colnames(df)) {
    stop("Column `", use, "` not found in ASE result.")
  }

  df[df[[use]] == geneID, , drop = FALSE]
}
