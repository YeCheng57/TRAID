#' Differential splicing analysis at event level
#'
#' Perform differential splicing analysis between two groups of samples
#' using event-level PSI (or delta PSI) values from an `ASResult` object.
#' Statistical testing is performed per event using either Wilcoxon rank-sum test
#' (recommended) or Student's t-test.
#'
#' @param x An `ASResult` object returned by `runAS()`.
#'
#' @param sample_info A data.frame containing sample metadata.
#' Must include a column `sample_id` and the grouping variable specified by `group_col`.
#'
#' @param group_col Character. Column name in `sample_info` defining the grouping variable.
#'
#' @param group1 Character. Name of the first group.
#'
#' @param group2 Character. Name of the second group.
#'
#' @param site Character. One of `"start"`, `"end"`, or `"both"`.
#' Determines which AS events are used:
#' \describe{
#'   \item{"start"}{Donor (5') events (`result_start`)}
#'   \item{"end"}{Acceptor (3') events (`result_end`)}
#'   \item{"both"}{Combined events from both}
#' }
#'
#' @param value Character. One of `"psi"` or `"delta_psi"`.
#' The metric used for comparison between groups.
#'
#' @param test Character. One of `"wilcox"` or `"t"`.
#' Statistical test used to compare groups.
#' Default is `"wilcox"` (recommended for PSI data).
#'
#' @param min_samples Integer. Minimum number of non-missing samples required
#' in each group for a given event. Default is 3.
#'
#' @return A data.frame containing differential splicing results at event level:
#' \describe{
#'   \item{event_id}{Splicing event identifier (chrom:start-end)}
#'   \item{chrom, start, end}{Genomic coordinates}
#'   \item{gene_name, gene_id}{Gene annotation (if available)}
#'   \item{event_source}{"start" or "end" (if `site = "both"`)}
#'   \item{mean_group1, mean_group2}{Mean PSI (or delta PSI) per group}
#'   \item{delta}{Difference between groups (group1 - group2)}
#'   \item{pvalue}{Raw p-value}
#'   \item{padj}{FDR-adjusted p-value}
#'   \item{direction}{Group with higher mean value}
#'   \item{n1, n2}{Number of samples per group}
#' }
#'
#' Results are sorted by absolute delta in descending order.
#'
#' @details
#' PSI values are bounded between 0 and 1 and often non-normally distributed.
#' Therefore, Wilcoxon rank-sum test is used by default.
#'
#' This function operates at the event level and does not perform gene-level aggregation.
#'
#' @export
findASDifference <- function(
    x,
    sample_info,
    group_col,
    group1,
    group2,
    site = c("start", "end", "both"),
    value = c("psi", "delta_psi"),
    test = c("wilcox", "t"),
    min_samples = 3
) {

  # ---- checks ----
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  if (!"sample_id" %in% colnames(sample_info)) {
    stop("`sample_info` must contain column `sample_id`.")
  }

  if (!group_col %in% colnames(sample_info)) {
    stop("`group_col` not found in sample_info.")
  }

  site <- match.arg(site)
  value <- match.arg(value)
  test <- match.arg(test)

  # ---- select data ----
  df <- switch(
    site,
    start = transform(x$result_start, event_source = "start"),
    end = transform(x$result_end, event_source = "end"),
    both = rbind(
      transform(x$result_start, event_source = "start"),
      transform(x$result_end, event_source = "end")
    )
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
  if (!"event_id" %in% colnames(df)) {
    df$event_id <- paste0(df$chrom, ":", df$start, "-", df$end)
  }

  # ---- merge metadata ----
  df <- merge(df, sample_info[, c("sample_id", group_col)], by = "sample_id")

  df <- df[df[[group_col]] %in% c(group1, group2), , drop = FALSE]

  if (nrow(df) == 0) {
    stop("No samples found for selected groups.")
  }

  # ---- split by event ----
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

    pval <- tryCatch({
      if (test == "wilcox") {
        stats::wilcox.test(g1, g2)$p.value
      } else {
        stats::t.test(g1, g2)$p.value
      }
    }, error = function(e) NA)

    data.frame(
      event_id = unique(d$event_id),

      chrom = unique(d$chrom),
      start = unique(d$start),
      end = unique(d$end),

      gene_name = if ("gene_names" %in% colnames(d))
        unique(d$gene_names)[1] else NA,

      gene_id = if ("gene_ids" %in% colnames(d))
        unique(d$gene_ids)[1] else NA,

      event_source = unique(d$event_source),

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

  # ---- multiple testing ----
  res$padj <- stats::p.adjust(res$pvalue, method = "fdr")

  # ---- direction ----
  res$direction <- ifelse(res$delta > 0, group1, group2)

  # ---- sort ----
  res <- res[order(abs(res$delta), decreasing = TRUE), ]

  rownames(res) <- NULL
  res
}
