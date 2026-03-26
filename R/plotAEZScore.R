#' Plot top AE z-scores for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Sample ID.
#' @param topN Number of top genes ranked by absolute z-score to display.
#' @param groupBy Optional grouping column in `sample_info`; only used for subtitle.
#' @param sample_info Optional sample metadata overriding `x$metadata$sample_info`.
#'
#' @return A ggplot object.
#' @export
plotAEZScore <- function(
    x,
    sampleID,
    topN = 30,
    groupBy = NULL,
    sample_info = NULL
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

  if (!is.null(topN) && is.finite(topN)) {
    topN <- as.integer(topN)
    if (topN <= 0) {
      stop("`topN` must be a positive integer.")
    }
    df <- utils::head(df, topN)
  }

  label_col <- if ("gene_name" %in% colnames(df)) "gene_name" else "gene_id"
  if (!label_col %in% colnames(df)) {
    df$label <- seq_len(nrow(df))
  } else {
    df$label <- df[[label_col]]
  }

  df <- df[order(df$zscore), , drop = FALSE]
  df$label <- factor(df$label, levels = df$label)

  si <- if (!is.null(sample_info)) sample_info else x$metadata$sample_info
  sample_group <- NULL
  if (!is.null(groupBy) && !is.null(si) && groupBy %in% colnames(si)) {
    m <- match(sampleID, si$sample_id)
    sample_group <- si[[groupBy]][m]
  }

  ggplot2::ggplot(df, ggplot2::aes(x = label, y = zscore)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("AE z-scores: ", sampleID),
      subtitle = if (!is.null(sample_group)) paste(groupBy, "=", sample_group) else NULL,
      x = NULL,
      y = "Z-score"
    )
}
