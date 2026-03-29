#' Summarize AEResult object
#'
#' @param object An `AEResult` object.
#' @param ... Reserved.
#'
#' @return A summary object.
#' @export
summary.AEResult <- function(object, ...) {
  if (!inherits(object, "AEResult")) {
    stop("`object` must be an AEResult.")
  }

  res <- object$result

  out <- list(
    n_samples = length(unique(res$sample_id)),
    n_genes = length(unique(res$gene_id)),
    n_rows = nrow(res),
    n_outliers = sum(res$is_outlier, na.rm = TRUE)
  )

  class(out) <- "summary.AEResult"
  out
}
#' @export
print.summary.AEResult <- function(x, ...) {
  cat("AEResult summary\n")
  cat(" Samples :", x$n_samples, "\n")
  cat(" Genes   :", x$n_genes, "\n")
  cat(" Rows    :", x$n_rows, "\n")
  cat(" Outliers:", x$n_outliers, "\n")
  invisible(x)
}

#' @export
plot.AEResult <- function(
    x,
    type = c("gene", "gene_oe", "sample_summary", "volcano", "zscore"),
    geneID = NULL,
    sampleID = NULL,
    groupBy = NULL,
    sample_info = NULL,
    ...
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  type <- match.arg(type)

  if (type == "gene") {
    if (is.null(geneID)) stop("`geneID` is required for `type = 'gene'`.")
    return(plotAEGene(x, geneID = geneID, groupBy = groupBy, sample_info = sample_info, ...))
  }

  if (type == "gene_oe") {
    if (is.null(geneID)) stop("`geneID` is required for `type = 'gene_oe'`.")
    return(plotAEGeneObservedExpected(x, geneID = geneID, groupBy = groupBy, sample_info = sample_info, ...))
  }

  if (type == "sample_summary") {
    return(plotAESampleSummary(x, groupBy = groupBy, sample_info = sample_info, ...))
  }

  if (type == "volcano") {
    if (is.null(sampleID)) stop("`sampleID` is required for `type = 'volcano'`.")
    return(plotAEVolcano(x, sampleID = sampleID, ...))
  }

  if (type == "zscore") {
    if (is.null(sampleID)) stop("`sampleID` is required for `type = 'zscore'`.")
    return(plotAEZScore(x, sampleID = sampleID, ...))
  }

  stop("Unsupported `type`: ", type)
}
