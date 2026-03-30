#' Find differential splicing events between two groups
#'
#' @param x An `ASResult` object.
#' @param sample_info Data.frame containing `sample_id` and grouping column.
#' @param group_col Column name in sample_info defining groups.
#' @param group1 First group.
#' @param group2 Second group.
#' @param site One of "start", "end", "both".
#' @param value "psi" or "delta_psi".
#' @param min_samples Minimum samples per group.
#'
#' @return Data.frame of differential splicing events.
#' @export
findASDifference <- function(
    x,
    sample_info,
    group_col,
    group1,
    group2,
    site = c("start", "end", "both"),
    value = c("psi", "delta_psi"),
    min_samples = 3
) {

  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  value <- match.arg(value)

  if (!"sample_id" %in% colnames(sample_info)) {
    stop("sample_info must contain `sample_id`.")
  }
  if (!group_col %in% colnames(sample_info)) {
    stop("group_col not found in sample_info.")
  }

  # ---- 取数据 ----
  df <- switch(
    site,
    start = x$result_start,
    end = x$result_end,
    both = rbind(x$result_start, x$result_end)
  )

  if (is.null(df) || nrow(df) == 0) {
    stop("No AS data available.")
  }

  req <- c("sample_id", "chrom", "start", "end", value)
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  # ---- event_id ----
  df$event_id <- paste0(df$chrom, ":", df$start, "-", df$end)

  # ---- merge 分组 ----
  df <- merge(df, sample_info[, c("sample_id", group_col)],
              by = "sample_id")

  # ---- 只保留两组 ----
  df <- df[df[[group_col]] %in% c(group1, group2), ]

  if (nrow(df) == 0) {
    stop("No samples in selected groups.")
  }

  # ---- 主循环（event-level）----
  split_list <- split(df, df$event_id)

  res_list <- lapply(split_list, function(d) {

    g1 <- d[d[[group_col]] == group1, value]
    g2 <- d[d[[group_col]] == group2, value]

    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]

    if (length(g1) < min_samples || length(g2) < min_samples) {
      return(NULL)
    }

    mean1 <- mean(g1)
    mean2 <- mean(g2)
    delta <- mean1 - mean2

    # Wilcoxon（对PSI更稳）
    pval <- tryCatch(
      stats::wilcox.test(g1, g2)$p.value,
      error = function(e) NA
    )

    data.frame(
      event_id = unique(d$event_id),
      mean_group1 = mean1,
      mean_group2 = mean2,
      delta = delta,
      pvalue = pval,
      n1 = length(g1),
      n2 = length(g2),
      stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, res_list)

  if (is.null(res) || nrow(res) == 0) {
    stop("No valid events found.")
  }

  # ---- 多重校正 ----
  res$padj <- stats::p.adjust(res$pvalue, method = "fdr")

  # ---- 排序 ----
  res <- res[order(abs(res$delta), decreasing = TRUE), ]

  rownames(res) <- NULL
  res
}
