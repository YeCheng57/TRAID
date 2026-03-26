#' Plot AE values for one gene across samples
#'
#' @param x An `AEResult` object.
#' @param geneID Gene identifier or gene symbol.
#' @param use One of `"gene_id"` or `"gene_name"`.
#' @param value Which value to plot. One of
#'   `"raw_count"`, `"norm_count"`, `"expected_count"`, `"zscore"`, `"log2fc"`.
#' @param logScale Logical; whether to use log10 scale for count-like values.
#' @param groupBy Optional grouping column in `sample_info`.
#' @param sample_info Optional sample metadata overriding `x$metadata$sample_info`.
#'
#' @return A ggplot object.
#' @export
plotAEGene <- function(
    x,
    geneID,
    use = c("gene_id", "gene_name"),
    value = c("raw_count", "norm_count", "expected_count", "zscore", "log2fc"),
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

  if (!use %in% colnames(df)) {
    stop("Column `", use, "` not found in AE result.")
  }
  if (!value %in% colnames(df)) {
    stop("Column `", value, "` not found in AE result.")
  }

  df <- df[df[[use]] == geneID, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No AE result found for gene: ", geneID)
  }

  df <- .attach_group_info(df, x, groupBy = groupBy, sample_info = sample_info)

  if (!"group" %in% colnames(df)) {
    df$group <- "All"
  }

  df <- df[order(df$group, df[[value]]), , drop = FALSE]
  df$sample_id <- factor(df$sample_id, levels = df$sample_id)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = sample_id, y = .data[[value]], fill = group)
  ) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("AE: ", geneID),
      x = "Sample",
      y = value
    )

  if (identical(value, "zscore")) {
    p <- p + ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed")
  }

  if (isTRUE(logScale)) {
    if (!value %in% c("raw_count", "norm_count", "expected_count")) {
      warning("`logScale = TRUE` is mainly intended for count-like values.")
    }
    p <- p + ggplot2::scale_y_continuous(trans = "log10")
  }

  p
}
