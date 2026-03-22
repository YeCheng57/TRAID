#' Plot a cartoon junction diagram for one sample versus controls
#'
#' @param x An `ASResult` object, preferably generated with `returnAll = TRUE`.
#' @param sampleID Character scalar specifying the target sample.
#' @param chrom Chromosome of the focal junction.
#' @param start Start coordinate of the focal junction.
#' @param end End coordinate of the focal junction.
#' @param site One of "start" or "end".
#' @param topN Integer; maximum number of junctions to show.
#'
#' @return A ggplot object.
#' @export
plotASJunctionCartoon <- function(
    x,
    sampleID,
    chrom,
    start,
    end,
    site = c("start", "end"),
    topN = 3
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Requested site result is empty.")
  }
  if (!sampleID %in% df$sample_id) {
    stop("`sampleID` not found in result.")
  }

  # 取同一事件簇
  if (site == "start") {
    event_df <- df[df$chrom == chrom & df$start == start, , drop = FALSE]
    partner_col <- "end"
    focal_partner <- end
    anchor_label <- paste0(chrom, ":", start)
  } else {
    event_df <- df[df$chrom == chrom & df$end == end, , drop = FALSE]
    partner_col <- "start"
    focal_partner <- start
    anchor_label <- paste0(chrom, ":", end)
  }

  if (nrow(event_df) == 0) {
    stop("Requested focal site not found.")
  }

  sample_df <- event_df[event_df$sample_id == sampleID, , drop = FALSE]
  others_df <- event_df[event_df$sample_id != sampleID, , drop = FALSE]

  if (nrow(sample_df) == 0) {
    stop("Target sample has no junctions at this focal site.")
  }

  # 其余样本统计
  if (nrow(others_df) > 0) {
    mean_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = others_df[[partner_col]]),
      FUN = mean,
      na.rm = TRUE
    )
    sd_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = others_df[[partner_col]]),
      FUN = stats::sd,
      na.rm = TRUE
    )
    colnames(mean_df)[2] <- "mean_count"
    colnames(sd_df)[2] <- "sd_count"
    others_sum <- merge(mean_df, sd_df, by = "partner", all = TRUE)
  } else {
    others_sum <- data.frame(
      partner = numeric(0),
      mean_count = numeric(0),
      sd_count = numeric(0)
    )
  }

  # 选展示 junction：focal 一定保留 + sample 中 count 最大若干条
  sample_rank <- sample_df[order(sample_df$junction_count, decreasing = TRUE), , drop = FALSE]
  keep_partners <- unique(c(focal_partner, head(sample_rank[[partner_col]], topN)))

  sample_df <- sample_df[sample_df[[partner_col]] %in% keep_partners, , drop = FALSE]
  others_sum <- others_sum[others_sum$partner %in% keep_partners, , drop = FALSE]

  # 固定全局 x 坐标 —— 上下两排完全共用
  partners <- sort(unique(keep_partners))
  x_anchor <- 0
  x_partner <- seq_along(partners) * 2.2 + 2.2
  names(x_partner) <- as.character(partners)

  # 行高度
  y_sample <- 2.0
  y_control <- 0.0

  # exon 尺寸
  exon_w <- 0.55
  exon_h <- 0.18

  # exon 数据：上下两排共用同一套 x
  exon_one_row <- data.frame(
    xmid = c(x_anchor, x_partner),
    exon_type = c("anchor", ifelse(partners == focal_partner, "focal", "other")),
    stringsAsFactors = FALSE
  )

  exon_sample <- exon_one_row
  exon_sample$panel <- "sample"
  exon_sample$y <- y_sample

  exon_control <- exon_one_row
  exon_control$panel <- "control"
  exon_control$y <- y_control

  exon_df <- rbind(exon_sample, exon_control)
  exon_df$xmin <- exon_df$xmid - exon_w
  exon_df$xmax <- exon_df$xmid + exon_w
  exon_df$ymin <- exon_df$y - exon_h / 2
  exon_df$ymax <- exon_df$y + exon_h / 2

  # 样本弧线
  sample_plot <- data.frame(
    x = x_anchor + exon_w * 0.9,
    xend = x_partner[as.character(sample_df[[partner_col]])] - exon_w * 0.9,
    y = y_sample + exon_h / 2,
    yend = y_sample + exon_h / 2,
    label = as.character(sample_df$junction_count),
    focal = sample_df[[partner_col]] == focal_partner,
    label_x = (x_anchor + x_partner[as.character(sample_df[[partner_col]])]) / 2,
    label_y = y_sample + 0.42,
    stringsAsFactors = FALSE
  )

  # 对照弧线
  others_plot <- data.frame(
    x = x_anchor + exon_w * 0.9,
    xend = x_partner[as.character(others_sum$partner)] - exon_w * 0.9,
    y = y_control + exon_h / 2,
    yend = y_control + exon_h / 2,
    label = paste0(round(others_sum$mean_count, 0), " ± ", round(others_sum$sd_count, 0)),
    focal = others_sum$partner == focal_partner,
    label_x = (x_anchor + x_partner[as.character(others_sum$partner)]) / 2,
    label_y = y_control + 0.42,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = exon_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = exon_type),
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_curve(
      data = sample_plot,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = focal),
      curvature = -0.25,
      linewidth = 0.6
    ) +
    ggplot2::geom_curve(
      data = others_plot,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = focal),
      curvature = -0.25,
      linewidth = 0.6
    ) +
    ggplot2::geom_text(
      data = sample_plot,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      size = 5
    ) +
    ggplot2::geom_text(
      data = others_plot,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      size = 5
    ) +
    ggplot2::annotate("text", x = x_anchor - 1.2, y = y_sample + 0.15, label = sampleID, hjust = 1, size = 6) +
    ggplot2::annotate("text", x = x_anchor - 1.2, y = y_control + 0.15, label = "Controls (mean ± SD)", hjust = 1, size = 6) +
    ggplot2::scale_fill_manual(
      values = c(anchor = "#2B1E1A", focal = "#FFFFFF", other = "#BDBDBD"),
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#7A5230", `FALSE` = "#7A5230"),
      guide = "none"
    ) +
    ggplot2::annotate("text", x = x_anchor, y = y_sample + 0.62, label = anchor_label, size = 4) +
    ggplot2::coord_cartesian(
      xlim = c(x_anchor - 1.8, max(x_partner) + 1.2),
      ylim = c(-0.5, 2.9),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(20, 30, 20, 90)
    )
}
