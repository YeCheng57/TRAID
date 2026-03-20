#' Plot expression pattern of one gene across samples
#'
#' @param x An `AEResult` object.
#' @param geneID Character scalar specifying the target gene.
#' @param sampleID Optional character scalar specifying one sample to highlight.
#' @param value Which value to plot. One of "raw_count", "norm_count", "expected_count", "zscore", "log2fc".
#' @param logScale Logical; whether to use log10 scale for count-based values.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#'
#' @return A ggplot object.
#' @export
plotAEGene <- function(
    x,
    geneID,
    sampleID = NULL,
    value = c("raw_count", "norm_count", "expected_count", "zscore", "log2fc"),
    logScale = TRUE,
    groupBy = NULL
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  value <- match.arg(value)
  df <- x$result

  if (!"gene_id" %in% colnames(df)) {
    stop("No `gene_id` column found in result.")
  }
  if (!geneID %in% df$gene_id) {
    stop("`geneID` not found in result.")
  }
  if (!value %in% colnames(df)) {
    stop("Requested `value` column not found in result.")
  }

  df <- df[df$gene_id == geneID, , drop = FALSE]
  df <- .attach_group_info(df, x, groupBy = groupBy)

  if (!is.null(sampleID)) {
    if (!sampleID %in% df$sample_id) {
      stop("`sampleID` not found for this gene.")
    }
    df$highlight <- df$sample_id == sampleID
  } else {
    df$highlight <- FALSE
  }

  gene_label <- geneID
  if ("gene_name" %in% colnames(df)) {
    gn <- unique(df$gene_name)
    gn <- gn[!is.na(gn) & gn != ""]
    if (length(gn) > 0) {
      gene_label <- paste0(gn[1], " (", geneID, ")")
    }
  }

  df <- .order_samples_for_plot(
    df,
    sample_col = "sample_id",
    value_col = value,
    group_col = if (!is.null(groupBy)) "group" else NULL,
    decreasing = TRUE
  )

  base_mapping <- ggplot2::aes(
    x = sample_id,
    y = .data[[value]]
  )

  if (!is.null(groupBy) && "group" %in% colnames(df) && any(!is.na(df$group))) {
    p <- ggplot2::ggplot(df, base_mapping) +
      ggplot2::geom_point(
        ggplot2::aes(color = group, shape = is_outlier),
        alpha = 0.7
      ) +
      ggplot2::labs(color = groupBy)
  } else {
    p <- ggplot2::ggplot(df, base_mapping) +
      ggplot2::geom_point(
        ggplot2::aes(color = is_outlier),
        alpha = 0.7
      )
  }

  if (!is.null(sampleID)) {
    p <- p +
      ggplot2::geom_point(
        data = df[df$highlight, , drop = FALSE],
        size = 3
      )
  }

  if (isTRUE(logScale) && value %in% c("raw_count", "norm_count", "expected_count")) {
    p <- p + ggplot2::scale_y_log10()
  }

  p +
    ggplot2::labs(
      title = paste("AE gene plot:", gene_label),
      x = "Sample",
      y = value
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}
