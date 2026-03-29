#' Plot observed and expected expression for one gene
#'
#' @param x An `AEResult` object.
#' @param geneID Gene identifier or gene symbol.
#' @param use One of `"gene_id"` or `"gene_name"`.
#' @param logScale Logical; whether to log10-scale the y axis.
#' @param groupBy Optional grouping column in `sample_info`.
#' @param sample_info Optional sample metadata overriding `x$metadata$sample_info`.
#' @param orderBy One of `"observed"`, `"expected"`, or `"difference"`.
#' @param decreasing Logical; whether to sort in decreasing order.
#' @param showSampleLabels Logical; whether to show sample labels on x axis.
#'
#' @return A ggplot object.
#' @export
plotAEGeneObservedExpected <- function(
    x,
    geneID,
    use = c("gene_id", "gene_name"),
    logScale = FALSE,
    groupBy = NULL,
    sample_info = NULL,
    orderBy = c("observed", "expected", "difference"),
    decreasing = TRUE,
    showSampleLabels = FALSE
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  use <- match.arg(use)
  value <- 'raw_count'
  orderBy <- match.arg(orderBy)

  df <- x$result

  if (!use %in% colnames(df)) {
    stop("Column `", use, "` not found in AE result.")
  }

  req <- c("sample_id", value, "expected_count")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns in AE result: ", paste(miss, collapse = ", "))
  }

  df <- df[df[[use]] == geneID, , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No AE result found for gene: ", geneID)
  }

  df <- .attach_group_info(df, x, groupBy = groupBy, sample_info = sample_info)

  if (orderBy == "observed") {
    ord_val <- df[[value]]
  } else if (orderBy == "expected") {
    ord_val <- df$expected_count
  } else {
    ord_val <- abs(df[[value]] - df$expected_count)
  }

  if ("group" %in% colnames(df)) {
    df <- df[order(df$group, ord_val, decreasing = decreasing), , drop = FALSE]
  } else {
    df <- df[order(ord_val, decreasing = decreasing), , drop = FALSE]
  }

  df$sample_id <- factor(df$sample_id, levels = df$sample_id)

  obs_df <- data.frame(
    sample_id = df$sample_id,
    value = df[[value]],
    type = "Observed",
    stringsAsFactors = FALSE
  )

  exp_df <- data.frame(
    sample_id = df$sample_id,
    value = df$expected_count,
    type = "Expected",
    stringsAsFactors = FALSE
  )

  if ("group" %in% colnames(df)) {
    obs_df$group <- df$group
    exp_df$group <- df$group
  }

  plot_df <- rbind(obs_df, exp_df)
  plot_df$type <- factor(plot_df$type, levels = c("Expected", "Observed"))

  if ("group" %in% colnames(plot_df)) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = sample_id, y = value, color = group, shape = type)
    ) +
      ggplot2::geom_point(size = 2.5)
  } else {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = sample_id, y = value, color = type)
    ) +
      ggplot2::geom_point(size = 2.5)
  }

  p <- p +
    ggplot2::scale_shape_manual(values = c(16, 17)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = if (showSampleLabels) {
        ggplot2::element_text(angle = 90, hjust = 1, size = 7)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (showSampleLabels) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      }
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
