#' Plot a sashimi-like junction view for one sample versus the remaining cohort
#'
#' @param x An `ASResult` object, preferably generated with `returnAll = TRUE`.
#' @param sampleID Character scalar specifying the target sample.
#' @param chrom Chromosome of the focal junction.
#' @param start Start coordinate of the focal junction.
#' @param end End coordinate of the focal junction.
#' @param site One of "start" or "end".
#' @param value One of "junction_count" or "psi".
#'
#' @return A ggplot object.
#' @export
plotASSashimi <- function(
    x,
    sampleID,
    chrom,
    start,
    end,
    site = c("start", "end"),
    value = c("junction_count", "psi")
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  value <- match.arg(value)

  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Requested site result is empty.")
  }

  if (!sampleID %in% df$sample_id) {
    stop("`sampleID` not found in result.")
  }

  req <- c("sample_id", "chrom", "start", "end", value)
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  # 确定 focal site
  if (site == "start") {
    focal_df <- df[df$chrom == chrom & df$start == start, , drop = FALSE]
    focal_pos <- start
    partner_col <- "end"
  } else {
    focal_df <- df[df$chrom == chrom & df$end == end, , drop = FALSE]
    focal_pos <- end
    partner_col <- "start"
  }

  if (nrow(focal_df) == 0) {
    stop("No junctions found for the requested focal site.")
  }

  # 只保留和 focal junction 同一个事件簇的连接
  # 若 site == start，则固定 chrom+start；
  # 若 site == end，则固定 chrom+end。
  event_df <- focal_df

  # 目标样本
  sample_df <- event_df[event_df$sample_id == sampleID, , drop = FALSE]
  if (nrow(sample_df) == 0) {
    stop("Target sample has no junctions at the requested focal site.")
  }

  # 其余样本
  others_df <- event_df[event_df$sample_id != sampleID, , drop = FALSE]

  # 对其余样本按 partner site 聚合 mean/sd
  if (nrow(others_df) > 0) {
    mean_df <- stats::aggregate(
      others_df[[value]],
      by = list(partner = others_df[[partner_col]]),
      FUN = mean,
      na.rm = TRUE
    )
    sd_df <- stats::aggregate(
      others_df[[value]],
      by = list(partner = others_df[[partner_col]]),
      FUN = stats::sd,
      na.rm = TRUE
    )

    colnames(mean_df)[2] <- "mean_value"
    colnames(sd_df)[2] <- "sd_value"

    others_sum <- merge(mean_df, sd_df, by = "partner", all = TRUE, sort = FALSE)
  } else {
    others_sum <- data.frame(
      partner = numeric(0),
      mean_value = numeric(0),
      sd_value = numeric(0)
    )
  }

  # 目标样本整理
  sample_plot <- data.frame(
    panel = "Sample",
    partner = sample_df[[partner_col]],
    value = sample_df[[value]],
    lower = NA_real_,
    upper = NA_real_,
    stringsAsFactors = FALSE
  )

  # 其余样本整理
  others_plot <- data.frame(
    panel = "Others",
    partner = others_sum$partner,
    value = others_sum$mean_value,
    lower = others_sum$mean_value - others_sum$sd_value,
    upper = others_sum$mean_value + others_sum$sd_value,
    stringsAsFactors = FALSE
  )

  plot_df <- rbind(sample_plot, others_plot)

  # 为了让图从左到右好看，按 partner 排序
  plot_df <- plot_df[order(plot_df$partner), , drop = FALSE]
  plot_df$panel <- factor(plot_df$panel, levels = c("Sample", "Others"))

  # 生成弧线高度：按 value 映射
  # 用 curve 表示连接，点表示连接终点，竖线表示其他样本 mean ± sd
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = focal_pos, xend = partner, y = 0, yend = value)) +
    ggplot2::geom_curve(
      curvature = 0.2,
      alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = partner, y = value),
      size = 2
    ) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = paste0("AS sashimi-like plot: ", chrom, ":", start, "-", end, " (", site, ")"),
      x = "Genomic position",
      y = value
    )

  # 只在 Others 面板加 mean ± sd 误差线
  if (nrow(others_plot) > 0) {
    p <- p +
      ggplot2::geom_errorbar(
        data = others_plot,
        ggplot2::aes(
          x = partner,
          ymin = lower,
          ymax = upper
        ),
        width = 0
      )
  }

  p
}
