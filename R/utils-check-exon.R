.check_exon_annotation <- function(x) {
  if (!is.data.frame(x)) {
    stop("`exon_annotation` must be a data.frame.")
  }

  required <- c(
    "chrom", "start", "end", "strand",
    "gene_id", "gene_name",
    "transcript_id", "exon_id"
  )

  miss <- setdiff(required, colnames(x))
  if (length(miss) > 0) {
    stop(
      "`exon_annotation` is missing required columns: ",
      paste(miss, collapse = ", ")
    )
  }

  if (!is.numeric(x$start)) {
    stop("`start` must be numeric/integer.")
  }
  if (!is.numeric(x$end)) {
    stop("`end` must be numeric/integer.")
  }

  invisible(TRUE)
}
