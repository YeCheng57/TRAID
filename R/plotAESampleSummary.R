#' Plot AE outlier burden per sample
#'
#' @param x An `AEResult` object.
#' @param groupBy Optional grouping column in `sample_info`.
#' @param sample_info Optional sample metadata overriding `x$metadata$sample_info`.
#' @param decreasing Logical; whether to sort in decreasing order of outlier burden.
#'
#' @return A ggplot object.
#' @export
plotAESampleSummary <- function(
    x,
    groupBy = NULL,
    sample_info = NULL,
    decreasing = TRUE
) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- x$result

  if (!all(c("sample_id", "is_outlier") %in% colnames(df))) {
    stop("AE result must contain `sample_id` and `is_outlier`.")
  }

  df_sum <- stats::aggregate(
    is_outlier ~ sample_id,
    data = df,
    FUN = function(z) sum(z, na.rm = TRUE)
  )
  colnames(df_sum)[2] <- "n_outlier"

  df_sum <- .attach_group_info(df_sum, x, groupBy = groupBy, sample_info = sample_info)

  if ("group" %in% colnames(df_sum)) {
    df_sum <- df_sum[order(df_sum$group, df_sum$n_outlier, decreasing = decreasing), , drop = FALSE]
  } else {
    df_sum <- df_sum[order(df_sum$n_outlier, decreasing = decreasing), , drop = FALSE]
  }

  df_sum$sample_id <- factor(df_sum$sample_id, levels = df_sum$sample_id)

  if ("group" %in% colnames(df_sum)) {
    p <- ggplot2::ggplot(df_sum, ggplot2::aes(x = sample_id, y = n_outlier, fill = group)) +
      ggplot2::geom_col()
  } else {
    p <- ggplot2::ggplot(df_sum, ggplot2::aes(x = sample_id, y = n_outlier)) +
      ggplot2::geom_col()
  }

  p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "AE outlier burden by sample",
      x = "Sample",
      y = "Number of outlier genes"
    )
}
