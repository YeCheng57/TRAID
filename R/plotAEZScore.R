#' Plot top AE z-scores for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Sample ID.
#' @param topN Number of top genes ranked by absolute z-score to display.
#'
#' @return A ggplot object.
#' @export
plotAEZScore <- function(
    x,
    sampleID,
    topN = 30
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- x$result
  df <- df[df$sample_id == sampleID, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No AE result found for sample: ", sampleID)
  }

  if (!"zscore" %in% colnames(df)) {
    stop("AE result must contain `zscore`.")
  }

  df$abs_zscore <- abs(df$zscore)
  df <- df[order(df$abs_zscore, decreasing = TRUE), , drop = FALSE]

  if (!is.null(topN)) {
    topN <- as.integer(topN)
    if (topN <= 0) {
      stop("`topN` must be a positive integer.")
    }
    df <- utils::head(df, topN)
  }

  label_col <- if ("gene_name" %in% colnames(df)) "gene_name" else "gene_id"
  if (!label_col %in% colnames(df)) {
    stop("Neither `gene_name` nor `gene_id` found in AE result.")
  }

  df$label <- df[[label_col]]
  df <- df[order(df$zscore), , drop = FALSE]
  df$label <- factor(df$label, levels = df$label)

  ggplot2::ggplot(df, ggplot2::aes(x = zscore, y = label)) +
    ggplot2::geom_col() +
    ggplot2::geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = paste0("Top AE z-scores: ", sampleID),
      x = "Z-score",
      y = NULL
    )
}
