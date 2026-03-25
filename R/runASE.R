#' Run site-level ASE analysis
#'
#' @param x A merged ASEReadCounter table.
#' @param annotation Optional transcript annotation list with components `exon` and `utr`.
#' @param gtf Optional local GTF path. Used only if `annotation` is not provided.
#' @param minDP Minimum total read depth.
#' @param region_keep Regions to retain, default `c("exonic", "UTR")`.
#' @param pCutoff Adjusted p-value cutoff.
#' @param ratioCutoff Allelic ratio cutoff. Sites with `alt_ratio > ratioCutoff`
#' or `< 1 - ratioCutoff` are considered significant.
#' @param returnAll Logical; whether to return all tested sites. If FALSE, only outliers are returned.
#'
#' @return An object of class `ASEResult`.
#' @export
runASE <- function(
    x,
    annotation = NULL,
    gtf = NULL,
    minDP = 10,
    region_keep = c("exonic", "UTR"),
    pCutoff = 0.05,
    ratioCutoff = 0.7,
    returnAll = FALSE
) {
  df <- parseASEReadCounter(x)

  # 只保留SNV
  df <- df[df$variant_type == "SNV", , drop = FALSE]

  # 注释 region
  if (!"region" %in% colnames(df)) {
    if (is.null(annotation)) {
      if (is.null(gtf)) {
        stop("Please provide either `annotation`, `gtf`, or a pre-annotated `region` column.")
      }
      annotation <- buildTranscriptAnnotation(gtf)
    }
    df <- .annotateASERegion(df, annotation)
  }

  # DP & region filter
  df <- df[df$total_count >= minDP, , drop = FALSE]
  df <- df[df$region %in% region_keep, , drop = FALSE]

  # 位点级统计
  df <- .run_binom_ase(df, pAdjustMethod = "fdr")

  # outlier
  df$is_outlier <- !is.na(df$padj) &
    df$padj < pCutoff &
    (df$alt_ratio > ratioCutoff | df$alt_ratio < (1 - ratioCutoff))

  # 保留列
  keep_cols <- c(
    "sample_id",
    "contig",
    "position",
    "ref",
    "alt",
    "variant_id",
    "variant_type",
    "ref_count",
    "alt_count",
    "total_count",
    "alt_ratio",
    "major_ratio",
    "region",
    "gene_name",
    "gene_id",
    "pvalue",
    "padj",
    "is_outlier"
  )
  keep_cols <- intersect(keep_cols, colnames(df))
  df <- df[, keep_cols, drop = FALSE]

  if (!isTRUE(returnAll)) {
    df <- df[!is.na(df$is_outlier) & df$is_outlier, , drop = FALSE]
  }

  rownames(df) <- NULL

  new_ASEResult(
    result = df,
    metadata = list(
      n_samples = length(unique(df$sample_id)),
      returnAll = returnAll
    ),
    parameters = list(
      minDP = minDP,
      region_keep = region_keep,
      pCutoff = pCutoff,
      ratioCutoff = ratioCutoff
    )
  )
}
