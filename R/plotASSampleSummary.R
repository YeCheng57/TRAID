#' Plot AS outlier counts per sample
#'
#' @param x An `ASResult` object.
#' @param site One of "start" or "end".
#' @param groupBy Optional character scalar specifying a column name in stored `colData`.
#'
#' @return A ggplot object.
#' @export
plotASSampleSummary <- function(x, site = c("start", "end"), groupBy = NULL) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Requested site result is empty.")
  }

  summ <- stats::aggregate(
    is_outlier ~ sample_id,
    data = df,
    FUN = function(z) sum(z, na.rm = TRUE)
  )
  colnames(summ)[2] <- "n_outlier"

  summ <- .attach_group_info(summ, x, groupBy = groupBy)

  if (!is.null(groupBy) && "group" %in% colnames(summ) && any(!is.na(summ$group))) {
    ggplot2::ggplot(
      summ,
      ggplot2::aes(
        x = stats::reorder(sample_id, -n_outlier),
        y = n_outlier,
        fill = group
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        title = paste("AS outlier counts per sample -", site),
        x = "Sample",
        y = "Number of outlier junctions",
        fill = groupBy
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  } else {
    ggplot2::ggplot(
      summ,
      ggplot2::aes(
        x = stats::reorder(sample_id, -n_outlier),
        y = n_outlier
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        title = paste("AS outlier counts per sample -", site),
        x = "Sample",
        y = "Number of outlier junctions"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
      )
  }
}
