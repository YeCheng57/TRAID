#' Plot AS z-score ranking for one sample
#'
#' @param x An `ASResult` object.
#' @param sampleID Character scalar specifying the target sample.
#' @param site One of "start" or "end".
#' @param topN Integer; number of top junctions to label.
#'
#' @return A ggplot object.
#' @export
plotASZScore <- function(x, sampleID, site = c("start", "end"), topN = 10) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  df <- getSampleAS(x, sampleID, site = site)

  if (nrow(df) == 0) {
    stop("No rows found for this sample.")
  }
  if (!"zscore" %in% colnames(df)) {
    stop("Result must contain `zscore` column.")
  }

  df <- df[order(abs(df$zscore), decreasing = TRUE), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df$label <- ""

  idx <- utils::head(seq_len(nrow(df)), topN)
  if (length(idx) > 0) {
    if ("gene_names" %in% colnames(df)) {
      lab <- df$gene_names[idx]
      bad <- is.na(lab) | lab == ""
      lab[bad] <- paste0(df$chrom[idx][bad], ":", df$start[idx][bad], "-", df$end[idx][bad])
      df$label[idx] <- lab
    } else {
      df$label[idx] <- paste0(df$chrom[idx], ":", df$start[idx], "-", df$end[idx])
    }
  }

  ggplot2::ggplot(df, ggplot2::aes(x = rank, y = zscore)) +
    ggplot2::geom_point(ggplot2::aes(color = is_outlier), alpha = 0.7) +
    ggplot2::geom_text(
      data = df[df$label != "", , drop = FALSE],
      ggplot2::aes(label = label),
      size = 3,
      vjust = -0.5,
      check_overlap = TRUE
    ) +
    ggplot2::labs(
      title = paste("AS z-score ranking:", sampleID, "-", site),
      x = "Rank by |z-score|",
      y = "z-score"
    )
}
