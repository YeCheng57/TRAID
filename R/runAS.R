#' Run aberrant splicing detection
#'
#' Perform junction-centric aberrant splicing analysis using regtools annotate results.
#'
#' @param x A regtools annotate table.
#' @param mode One of "both", "donor", "acceptor".
#' @param pAdjustMethod Multiple testing correction method.
#' @param pCutoff Adjusted p-value cutoff.
#' @param minTotalReads Minimum total reads at donor/acceptor site.
#' @param minJunctionReads Minimum reads supporting focal junction.
#'
#' @return An object of class `ASResult`.
#' @export
runAS <- function(
    x,
    mode = c("both", "donor", "acceptor"),
    pAdjustMethod = "bonferroni",
    pCutoff = 0.05,
    minTotalReads = 10,
    minJunctionReads = 5
) {
  mode <- match.arg(mode)

  df <- parseRegtoolsAnnotate(x)

  run_one_mode <- function(df, site) {
    out <- .compute_junction_psi(df, site = site)
    out <- .compute_as_statistics_leave_one_out(out, pAdjustMethod = pAdjustMethod)

    out$is_outlier <- !is.na(out$padj) &
      out$padj < pCutoff &
      out$total_site_count >= minTotalReads &
      out$junction_count >= minJunctionReads

    preferred <- c(
      "sample_id",
      "chrom",
      "start",
      "end",
      "strand",
      "splice_site",
      "anchor",
      "gene_names",
      "gene_ids",
      "transcripts",
      "junction_id",
      "site_type",
      "site_id",
      "junction_count",
      "total_site_count",
      "psi",
      "psi_mean_others",
      "psi_sd_others",
      "delta_psi",
      "zscore",
      "pvalue",
      "padj",
      "is_outlier"
    )

    keep_cols <- intersect(preferred, colnames(out))
    other_cols <- setdiff(colnames(out), keep_cols)
    out <- out[, c(keep_cols, other_cols), drop = FALSE]

    out
  }

  result_donor <- NULL
  result_acceptor <- NULL

  if (mode %in% c("both", "donor")) {
    result_donor <- run_one_mode(df, site = "donor")
  }

  if (mode %in% c("both", "acceptor")) {
    result_acceptor <- run_one_mode(df, site = "acceptor")
  }

  new_ASResult(
    result_donor = result_donor,
    result_acceptor = result_acceptor,
    metadata = list(
      mode = mode,
      n_rows_input = nrow(df),
      n_samples = length(unique(df$sample_id))
    ),
    parameters = list(
      pAdjustMethod = pAdjustMethod,
      pCutoff = pCutoff,
      minTotalReads = minTotalReads,
      minJunctionReads = minJunctionReads
    )
  )
}
