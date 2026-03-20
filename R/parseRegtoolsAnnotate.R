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
    "Sample", "chrom", "start", "end", "name", "score", "strand",
    "splice_site", "anchor", "gene_names", "gene_ids", "transcripts"
  )
  miss <- setdiff(req, colnames(x))
  if (length(miss) > 0) {
    stop("Missing required columns in regtools annotate table: ",
         paste(miss, collapse = ", "))
  }

  out <- data.frame(
    sample_id = as.character(x$Sample),
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

  out$junction_id <- paste0(out$chrom, ":", out$start, "-", out$end)
  out$donor_id <- paste(out$gene_names, out$chrom, out$start, sep = ":")
  out$acceptor_id <- paste(out$gene_names, out$chrom, out$end, sep = ":")

  out
}
