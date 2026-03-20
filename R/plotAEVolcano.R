#' Plot AE volcano plot for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Character scalar specifying the target sample.
#' @param topN Integer; number of top genes to label.
#'
#' @return A ggplot object.
#' @export
plotAEVolcano <- function(x, sampleID, topN = 10) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- getSampleAE(x, sampleID)

  if (!"padj" %in% colnames(df) || !"log2fc" %in% colnames(df)) {
    stop("Result must contain `padj` and `log2fc` columns.")
  }

  df$neglog10_padj <- -log10(df$padj)
  df$neglog10_padj[!is.finite(df$neglog10_padj)] <- NA_real_

  df$label <- ""
  ord <- order(df$padj, na.last = NA)
  if (length(ord) > 0) {
    idx <- utils::head(ord, topN)
    if ("gene_name" %in% colnames(df)) {
      lab <- df$gene_name[idx]
      bad <- is.na(lab) | lab == ""
      lab[bad] <- df$gene_id[idx][bad]
      df$label[idx] <- lab
    } else {
      df$label[idx] <- df$gene_id[idx]
    }
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = log2fc, y = neglog10_padj)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = is_outlier),
      alpha = 0.7
    ) +
    ggplot2::geom_text(
      data = df[df$label != "", , drop = FALSE],
      ggplot2::aes(label = label),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    ggplot2::labs(
      title = paste("AE volcano:", sampleID),
      x = "log2 fold change",
      y = expression(-log[10](adjusted~p))
    )
}
