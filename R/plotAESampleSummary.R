#' Plot AE outlier counts per sample
#'
#' @param x An `AEResult` object.
#'
#' @return A ggplot object.
#' @export
plotAESampleSummary <- function(x) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }

  df <- x$result

  if (!all(c("sample_id", "is_outlier") %in% colnames(df))) {
    stop("Result must contain `sample_id` and `is_outlier` columns.")
  }

  summ <- stats::aggregate(
    is_outlier ~ sample_id,
    data = df,
    FUN = function(z) sum(z, na.rm = TRUE)
  )
  colnames(summ)[2] <- "n_outlier"

  ggplot2::ggplot(
    summ,
    ggplot2::aes(x = stats::reorder(sample_id, -n_outlier), y = n_outlier)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = "AE outlier counts per sample",
      x = "Sample",
      y = "Number of outlier genes"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}
