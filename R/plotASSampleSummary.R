#' Plot AS outlier burden per sample
#'
#' Visualize the number or fraction of aberrant splicing events per sample.
#'
#' @param x An `ASResult` object.
#' @param site One of `"start"`, `"end"`, or `"both"`.
#' `"start"` uses `x$result_start`; `"end"` uses `x$result_end`;
#' `"both"` combines both tables.
#' @param groupBy Optional character; column name in `sample_info` used for grouping.
#' @param sample_info Optional data.frame; user-supplied sample metadata.
#' Must contain a column `sample_id`. If provided, it overrides
#' `x$metadata$sample_info`.
#' @param mode One of `"count"` or `"fraction"`.
#' @param decreasing Logical; whether to sort samples in decreasing order.
#'
#' @return A ggplot object.
#' @export
plotASSampleSummary <- function(
    x,
    site = c("start", "end", "both"),
    groupBy = NULL,
    sample_info = NULL,
    mode = c("count", "fraction"),
    decreasing = TRUE
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  mode <- match.arg(mode)

  if (site == "start") {
    df <- x$result_start
  } else if (site == "end") {
    df <- x$result_end
  } else {
    df <- rbind(
      x$result_start,
      x$result_end
    )
  }

  if (is.null(df) || nrow(df) == 0) {
    stop("No AS results available for the selected `site`.")
  }

  if (!all(c("sample_id", "is_outlier") %in% colnames(df))) {
    stop("AS result must contain `sample_id` and `is_outlier`.")
  }

  # 每个 sample 的 outlier 数
  df_count <- stats::aggregate(
    is_outlier ~ sample_id,
    data = df,
    FUN = function(z) sum(z, na.rm = TRUE)
  )
  colnames(df_count)[2] <- "n_outlier"

  # 每个 sample 的总事件数
  df_total <- stats::aggregate(
    is_outlier ~ sample_id,
    data = df,
    FUN = length
  )
  colnames(df_total)[2] <- "n_total"

  df_sum <- merge(df_count, df_total, by = "sample_id", all = TRUE)
  df_sum$fraction <- df_sum$n_outlier / df_sum$n_total

  # attach group
  df_sum <- .attach_group_info(
    df_sum,
    x,
    groupBy = groupBy,
    sample_info = sample_info
  )

  y_col <- if (mode == "count") "n_outlier" else "fraction"
  y_lab <- if (mode == "count") {
    "Number of aberrant splicing events"
  } else {
    "Fraction of aberrant splicing events"
  }

  # 排序
  if ("group" %in% colnames(df_sum)) {
    df_sum <- df_sum[order(df_sum$group, df_sum[[y_col]], decreasing = decreasing), , drop = FALSE]
  } else {
    df_sum <- df_sum[order(df_sum[[y_col]], decreasing = decreasing), , drop = FALSE]
  }

  df_sum$sample_id <- factor(df_sum$sample_id, levels = df_sum$sample_id)

  # 作图
  if ("group" %in% colnames(df_sum)) {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(x = sample_id, y = .data[[y_col]], fill = group)
    ) +
      ggplot2::geom_col()
  } else {
    p <- ggplot2::ggplot(
      df_sum,
      ggplot2::aes(x = sample_id, y = .data[[y_col]])
    ) +
      ggplot2::geom_col()
  }

  p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("AS outlier burden by sample (", site, ")"),
      x = "Sample",
      y = y_lab
    )
}
