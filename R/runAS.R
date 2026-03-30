#' Run aberrant splicing detection
#'
#' Perform junction-centric aberrant splicing analysis using regtools annotate results.
#'
#' @param x A regtools annotate table.
#' @param mode One of "both", "start", "end".
#' @param pAdjustMethod Multiple testing correction method.
#' @param pCutoff Adjusted p-value cutoff.
#' @param minTotalReads Minimum total reads at start/end site for defining outliers.
#' @param minJunctionReads Minimum reads supporting focal junction for defining outliers.
#' @param deltaPSI_cutoff Numeric; threshold for defining aberrant splicing events.
#' Default is 0.1 (10%). Events with absolute delta PSI ≥ this cutoff are considered aberrant.
#' @param sample_info Optional sample metadata as a data.frame.
#' Must contain a column `sample_id`, and can include additional annotations
#' such as `sex`, `Gene_diagnostic`, or user-defined grouping variables.
#' @param returnAll Logical; whether to return all tested junctions. If FALSE,
#' only outliers are returned.
#'
#' @return An object of class `ASResult`.
#' @export
runAS <- function(
    x,
    mode = c("both", "start", "end"),
    pAdjustMethod = "bonferroni",
    pCutoff = 0.05,
    minTotalReads = 10,
    minJunctionReads = 5,
    deltaPSI_cutoff = 0.1,
    sample_info = NULL,
    returnAll = FALSE
) {
  mode <- match.arg(mode)

  df <- parseRegtoolsAnnotate(x)
  sample_ids<-df$sample_id
  validateSamplesAgainstInfo(sample_ids, sample_info, input_name = "AE input")

  run_one_mode <- function(df, site) {
    # 先对所有 junction 计算 psi 和统计量
    out <- .compute_junction_psi(df, site = site)
    out <- .compute_as_statistics_leave_one_out(out, pAdjustMethod = pAdjustMethod)

    # 定义 outlier（阈值只影响 outlier 判定和默认输出，不影响 psi 计算）
    out$is_outlier <- !is.na(out$padj) &
      out$padj < pCutoff &
      out$total_site_count >= minTotalReads &
      out$junction_count >= minJunctionReads &
      out$delta_psi >= deltaPSI_cutoff

    # 如果统计上 NA，再补 reason_na
    low_total <- out$total_site_count < minTotalReads
    low_junc <- out$junction_count < minJunctionReads

    out$reason_na[is.na(out$reason_na) & low_total] <- "low_total_reads"
    out$reason_na[is.na(out$reason_na) & low_junc] <- "low_junction_reads"

    # 简化输出列
    keep_cols <- c(
      "sample_id",
      "chrom",
      "start",
      "end",
      "strand",
      "splice_site",
      "anchor",
      "gene_names",
      "transcripts",
      "site_type",
      "junction_count",
      "total_site_count",
      "psi",
      "psi_mean_others",
      "psi_sd_others",
      "delta_psi",
      "n_others_non_missing",
      "zscore",
      "pvalue",
      "padj",
      "reason_na",
      "is_outlier"
    )
    keep_cols <- intersect(keep_cols, colnames(out))
    out <- out[, keep_cols, drop = FALSE]

    # 默认只返回 outlier
    if (!isTRUE(returnAll)) {
      out <- out[!is.na(out$is_outlier) & out$is_outlier, , drop = FALSE]
    }

    rownames(out) <- NULL
    out
  }

  result_start <- NULL
  result_end <- NULL

  if (mode %in% c("both", "start")) {
    result_start <- run_one_mode(df, site = "start")
  }

  if (mode %in% c("both", "end")) {
    result_end <- run_one_mode(df, site = "end")
  }

  new_ASResult(
    result_start = result_start,
    result_end = result_end,
    metadata = list(
      mode = mode,
      sample_info = sample_info,
      n_rows_input = nrow(df),
      n_samples = length(unique(df$sample_id)),
      returnAll = returnAll
    ),
    parameters = list(
      pAdjustMethod = pAdjustMethod,
      pCutoff = pCutoff,
      minTotalReads = minTotalReads,
      minJunctionReads = minJunctionReads
    )
  )
}
