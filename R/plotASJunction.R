#' Plot PSI of one junction across samples
#'
#' @param x An `ASResult` object.
#' @param chrom Chromosome.
#' @param start Junction start.
#' @param end Junction end.
#' @param site One of "start" or "end".
#' @param sampleID Optional sample to highlight.
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#'
#' @return A ggplot object.
#' @export
plotASJunction <- function(
    x,
    chrom,
    start,
    end,
    site = c("start", "end"),
    sampleID = NULL,
    groupBy = NULL
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Requested site result is empty.")
  }

  df <- df[
    df$chrom == chrom &
      df$start == start &
      df$end == end,
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    stop("Requested junction not found.")
  }

  df <- .attach_group_info(df, x, groupBy = groupBy)

  if (!is.null(sampleID)) {
    if (!sampleID %in% df$sample_id) {
      stop("`sampleID` not found for this junction.")
    }
    df$highlight <- df$sample_id == sampleID
  } else {
    df$highlight <- FALSE
  }

  df <- .order_samples_for_plot(
    df,
    sample_col = "sample_id",
    value_col = "psi",
    group_col = if (!is.null(groupBy)) "group" else NULL,
    decreasing = TRUE
  )

  if (!is.null(groupBy) && "group" %in% colnames(df) && any(!is.na(df$group))) {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = sample_id, y = psi)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = group, shape = is_outlier),
        alpha = 0.8
      ) +
      ggplot2::labs(color = groupBy)
  } else {
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = sample_id, y = psi)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = is_outlier),
        alpha = 0.8
      )
  }

  if (!is.null(sampleID)) {
    p <- p +
      ggplot2::geom_point(
        data = df[df$highlight, , drop = FALSE],
        size = 3
      )
  }

  p +
    ggplot2::labs(
      title = paste0("AS junction PSI: ", chrom, ":", start, "-", end, " (", site, ")"),
      x = "Sample",
      y = "PSI"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}
