#' Plot observed vs expected expression for one gene
#'
#' @param x An `AEResult` object.
#' @param geneID Gene identifier or symbol.
#' @param use `"gene_id"` or `"gene_name"`.
#' @param value `"norm_count"` or `"raw_count"`
#' @param logScale Whether to log-scale y axis.
#' @param groupBy Grouping column.
#' @param sample_info Optional override metadata.
#'
#' @return ggplot object
#' @export
plotAEGeneObservedExpected <- function(
    x,
    geneID,
    use = c("gene_id", "gene_name"),
    value = c("norm_count", "raw_count"),
    logScale = FALSE,
    groupBy = NULL,
    sample_info = NULL
) {

  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  use <- match.arg(use)
  value <- match.arg(value)

  df <- x$result
  df <- df[df[[use]] == geneID, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No AE result found for gene: ", geneID)
  }

  df <- .attach_group_info(df, x, groupBy = groupBy, sample_info = sample_info)

  if (!"group" %in% colnames(df)) {
    df$group <- "All"
  }

  df <- df[order(df$group, df$sample_id), ]
  df$sample_id <- factor(df$sample_id, levels = df$sample_id)

  obs <- df[[value]]
  exp <- df$expected_count

  plot_df <- rbind(
    data.frame(sample_id = df$sample_id, group = df$group, value = obs, type = "Observed"),
    data.frame(sample_id = df$sample_id, group = df$group, value = exp, type = "Expected")
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = sample_id, y = value, color = group, shape = type)
  ) +
    ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("Observed vs Expected: ", geneID),
      x = "Sample",
      y = value
    )

  if (isTRUE(logScale)) {
    p <- p + ggplot2::scale_y_continuous(trans = "log10")
  }

  p
}
