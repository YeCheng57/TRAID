#' Plot observed vs expected expression for one gene across samples
#'
#' @param x An `AEResult` object.
#' @param geneID Character scalar specifying the target gene.
#' @param sampleID Optional character scalar specifying one sample to highlight.
#' @param logScale Logical; whether to use log10 scale.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#'
#' @return A ggplot object.
#' @export
plotAEGeneObservedExpected <- function(
    x,
    geneID,
    sampleID = NULL,
    logScale = TRUE,
    groupBy = NULL
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- x$result

  req <- c("gene_id", "sample_id", "raw_count", "expected_count")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  if (!geneID %in% df$gene_id) {
    stop("`geneID` not found in result.")
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
    value_col = "raw_count",
    group_col = if (!is.null(groupBy)) "group" else NULL,
    decreasing = TRUE
  )

  obs <- data.frame(
    sample_id = df$sample_id,
    value = df$raw_count,
    type = rep("Observed", nrow(df)),
    highlight = df$highlight,
    is_outlier = df$is_outlier,
    group = df$group,
    stringsAsFactors = FALSE
  )

  exp <- data.frame(
    sample_id = df$sample_id,
    value = df$expected_count,
    type = rep("Expected", nrow(df)),
    highlight = df$highlight,
    is_outlier = df$is_outlier,
    group = df$group,
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(obs, exp)
  plot_df$sample_id <- factor(plot_df$sample_id, levels = levels(df$sample_id))

  if (!is.null(groupBy) && any(!is.na(plot_df$group))) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = sample_id, y = value, shape = type, color = group)
    ) +
      ggplot2::geom_point(
        alpha = 0.8,
        position = ggplot2::position_dodge(width = 0.4)
      ) +
      ggplot2::labs(color = groupBy)
  } else {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = sample_id, y = value, shape = type, color = type)
    ) +
      ggplot2::geom_point(
        alpha = 0.8,
        position = ggplot2::position_dodge(width = 0.4)
      )
  }

  if (!is.null(sampleID)) {
    p <- p +
      ggplot2::geom_point(
        data = plot_df[plot_df$highlight, , drop = FALSE],
        size = 3,
        position = ggplot2::position_dodge(width = 0.4)
      )
  }

  if (isTRUE(logScale)) {
    p <- p + ggplot2::scale_y_log10()
  }

  p +
    ggplot2::labs(
      title = paste("Observed vs expected:", gene_label),
      x = "Sample",
      y = "Expression"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}
