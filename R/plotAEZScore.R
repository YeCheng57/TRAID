#' Plot AE z-score ranking for one sample
#'
#' @param x An `AEResult` object.
#' @param sampleID Character scalar specifying the target sample.
#' @param topN Integer; number of top genes to label.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#'
#' @return A ggplot object.
#' @export
plotAEZScore <- function(x, sampleID, topN = 10, groupBy = NULL) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- getSampleAE(x, sampleID)

  if (!"zscore" %in% colnames(df)) {
    stop("Result must contain `zscore` column.")
  }

  df <- df[order(abs(df$zscore), decreasing = TRUE), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df <- .attach_group_info(df, x, groupBy = groupBy)

  df$label <- ""
  idx <- utils::head(seq_len(nrow(df)), topN)
  if (length(idx) > 0) {
    if ("gene_name" %in% colnames(df)) {
      lab <- df$gene_name[idx]
      bad <- is.na(lab) | lab == ""
      lab[bad] <- df$gene_id[idx][bad]
      df$label[idx] <- lab
    } else {
      df$label[idx] <- df$gene_id[idx]
    }
  }

  if (!is.null(groupBy) && "group" %in% colnames(df) && any(!is.na(df$group))) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = rank, y = zscore)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = group, shape = is_outlier),
        alpha = 0.7
      )
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = rank, y = zscore)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = is_outlier),
        alpha = 0.7
      )
  }

  p +
    ggplot2::geom_text(
      data = df[df$label != "", , drop = FALSE],
      ggplot2::aes(label = label),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    ggplot2::labs(
      title = paste("AE z-score ranking:", sampleID),
      x = "Rank by |z-score|",
      y = "z-score"
    )
}
