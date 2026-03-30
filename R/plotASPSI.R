#' Plot PSI for a specific splice junction
#'
#' Visualize PSI values across samples for a given splice junction,
#' defined by chrom + start + end.
#'
#' @param x An `ASResult` object.
#'
#' @param chrom Chromosome.
#' @param start Junction start position.
#' @param end Junction end position.
#'
#' @param site One of `"start"` or `"end"`.
#' Determines whether to use `result_start` or `result_end`.
#'
#' @param gene Optional gene name to further filter events.
#'
#' @param groupBy Optional grouping column in `sample_info`.
#'
#' @param sample_info Optional data.frame; overrides `x$metadata$sample_info`.
#'
#' @param decreasing Logical; whether to sort samples by PSI.
#'
#' @return A ggplot object.
#' @export
plotASPSI <- function(
    x,
    chrom,
    start,
    end,
    site = c("start", "end"),
    gene = NULL,
    groupBy = NULL,
    sample_info = NULL,
    decreasing = TRUE
) {

  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)

  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Selected AS result table is empty.")
  }

  req <- c("sample_id", "chrom", "start", "end", "psi")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns in AS result: ",
         paste(miss, collapse = ", "))
  }

  # ---- 精确筛选 junction ----
  df <- df[
    df$chrom == chrom &
      df$start == start &
      df$end == end,
    ,
    drop = FALSE
  ]

  # gene 限制（可选）
  if (!is.null(gene) && "gene_names" %in% colnames(df)) {
    df <- df[df$gene_names == gene, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No AS junction found for specified coordinates.")
  }

  # 去 NA
  df <- df[!is.na(df$psi), , drop = FALSE]

  if (nrow(df) == 0) {
    stop("All PSI values are NA for this junction.")
  }

  # ---- attach group ----
  df <- .attach_group_info(
    df,
    x,
    groupBy = groupBy,
    sample_info = sample_info
  )

  # ---- 标题 ----
  event_label <- paste0(chrom, ":", start, "-", end)
  if (!is.null(gene)) {
    event_label <- paste0(gene, " | ", event_label)
  }

  # ---- 有分组：boxplot ----
  if ("group" %in% colnames(df)) {

    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = group, y = psi, fill = group)
    ) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.2, size = 2) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = paste0("AS PSI: ", event_label),
        x = groupBy,
        y = "PSI"
      )

  } else {

    # ---- 单样本排序 ----
    df <- df[order(df$psi, decreasing = decreasing), , drop = FALSE]
    df$sample_id <- factor(df$sample_id, levels = df$sample_id)

    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = sample_id, y = psi)
    ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = paste0("AS PSI: ", event_label),
        x = "Sample",
        y = "PSI"
      )
  }

  p
}
