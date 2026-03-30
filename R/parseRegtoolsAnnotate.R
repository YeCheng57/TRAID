#' Parse regtools annotate table
#'
#' @param x A regtools annotate result table.
#'
#' @return A standardized data.frame.
#' @keywords internal
parseRegtoolsAnnotate <- function(x) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.")
  }

  req <- c(
    "sample_id", "chrom", "start", "end", "name", "score", "strand",
    "splice_site", "anchor", "gene_names", "gene_ids", "transcripts"
  )
  miss <- setdiff(req, colnames(x))
  if (length(miss) > 0) {
    stop("Missing required columns in regtools annotate table: ",
         paste(miss, collapse = ", "))
  }

  out <- data.frame(
    sample_id = as.character(x$sample_id),
    chrom = as.character(x$chrom),
    start = as.integer(x$start),
    end = as.integer(x$end),
    junction_name = as.character(x$name),
    count = as.numeric(x$score),
    strand = as.character(x$strand),
    splice_site = as.character(x$splice_site),
    anchor = as.character(x$anchor),
    gene_names = as.character(x$gene_names),
    gene_ids = as.character(x$gene_ids),
    transcripts = as.character(x$transcripts),
    stringsAsFactors = FALSE
  )

  # ---- ⭐ 核心：collapse 相同 junction ----
  # key: sample_id + chrom + start + end + strand（建议保留 strand）

  out <- stats::aggregate(
    count ~ sample_id + chrom + start + end + strand,
    data = out,
    FUN = function(z) sum(z, na.rm = TRUE)
  )

  # ---- 恢复 annotation（取第一条即可）----
  meta_cols <- c("sample_id", "chrom", "start", "end", "strand")

  meta_df <- x[!duplicated(x[, meta_cols]), , drop = FALSE]

  meta_df <- data.frame(
    sample_id = as.character(meta_df$sample_id),
    chrom = as.character(meta_df$chrom),
    start = as.integer(meta_df$start),
    end = as.integer(meta_df$end),
    strand = as.character(meta_df$strand),
    junction_name = as.character(meta_df$name),
    splice_site = as.character(meta_df$splice_site),
    anchor = as.character(meta_df$anchor),
    gene_names = as.character(meta_df$gene_names),
    gene_ids = as.character(meta_df$gene_ids),
    transcripts = as.character(meta_df$transcripts),
    stringsAsFactors = FALSE
  )

  out <- merge(out, meta_df,
               by = c("sample_id", "chrom", "start", "end", "strand"),
               all.x = TRUE)

  # ---- 重新生成 ID ----
  out$junction_id <- paste0(out$chrom, ":", out$start, "-", out$end)

  out$start_id <- paste(out$gene_names, out$chrom, out$start, sep = ":")
  out$end_id <- paste(out$gene_names, out$chrom, out$end, sep = ":")

  out
}
