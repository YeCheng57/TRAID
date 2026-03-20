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

#' Plot AEResult object
#'
#' @param x An `AEResult` object.
#' @param type Type of plot: "volcano", "zscore", "sample_summary", "gene", "gene_oe".
#' @param sampleID Required for sample-level plots; optional for gene-level plots.
#' @param geneID Required for gene-level plots.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#' @param ... Additional arguments passed to underlying plotting functions.
#'
#' @return A ggplot object.
#' @export
plot.AEResult <- function(
    x,
    type = c("volcano", "zscore", "sample_summary", "gene", "gene_oe"),
    sampleID = NULL,
    geneID = NULL,
    groupBy = NULL,
    ...
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  type <- match.arg(type)

  if (type == "volcano") {
    if (is.null(sampleID)) stop("`sampleID` is required for volcano plot.")
    return(plotAEVolcano(x, sampleID = sampleID, groupBy = groupBy, ...))
  }

  if (type == "zscore") {
    if (is.null(sampleID)) stop("`sampleID` is required for z-score plot.")
    return(plotAEZScore(x, sampleID = sampleID, groupBy = groupBy, ...))
  }

  if (type == "sample_summary") {
    return(plotAESampleSummary(x, groupBy = groupBy, ...))
  }

  if (type == "gene") {
    if (is.null(geneID)) stop("`geneID` is required for gene plot.")
    return(plotAEGene(x, geneID = geneID, sampleID = sampleID, groupBy = groupBy, ...))
  }

  if (type == "gene_oe") {
    if (is.null(geneID)) stop("`geneID` is required for gene_oe plot.")
    return(plotAEGeneObservedExpected(x, geneID = geneID, sampleID = sampleID, groupBy = groupBy, ...))
  }

  stop("Unknown plot type.")
}
