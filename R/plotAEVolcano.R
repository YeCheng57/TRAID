#' Plot AE volcano plot for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Sample ID.
#' @param padjCutoff FDR cutoff for highlighting outliers.
#' @param zCutoff Optional z-score cutoff.
#' @param labelTopN Number of top genes (by |zscore|) to label.
#' @param groupBy Optional grouping column (only used in subtitle).
#' @param sample_info Optional sample metadata overriding stored metadata.
#'
#' @return ggplot object
#' @export
plotAEVolcano <- function(
    x,
    sampleID,
    padjCutoff = 0.05,
    zCutoff = NULL,
    labelTopN = 10,
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

  # ---- outlier定义 ----
  df$is_outlier <- !is.na(df$padj) & df$padj < padjCutoff
  if (!is.null(zCutoff)) {
    df$is_outlier <- df$is_outlier &
      !is.na(df$zscore) &
      abs(df$zscore) >= zCutoff
  }

  df$neglog10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))

  # ---- label ----
  df$abs_z <- abs(df$zscore)
  df <- df[order(df$abs_z, decreasing = TRUE), ]
  df$label <- ""

  if (!is.null(labelTopN) && labelTopN > 0) {
    top_idx <- seq_len(min(labelTopN, nrow(df)))
    label_col <- if ("gene_name" %in% colnames(df)) "gene_name" else "gene_id"
    df$label[top_idx] <- df[[label_col]][top_idx]
  }

  # ---- sample group（只用于subtitle）----
  si <- if (!is.null(sample_info)) sample_info else x$metadata$sample_info
  sample_group <- NULL

  if (!is.null(groupBy) && !is.null(si) && groupBy %in% colnames(si)) {
    m <- match(sampleID, si$sample_id)
    sample_group <- si[[groupBy]][m]
  }

  # ---- plot ----
  ggplot2::ggplot(df, ggplot2::aes(x = log2fc, y = neglog10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = is_outlier)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = Inf) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      title = paste0("AE Volcano: ", sampleID),
      subtitle = if (!is.null(sample_group)) paste(groupBy, "=", sample_group) else NULL,
      x = "log2 fold change",
      y = "-log10(FDR)"
    )
}
