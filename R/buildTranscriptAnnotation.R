#' Build transcript annotation from GTF
#'
#' @param gtf Path to a local GTF/GTF.GZ file.
#'
#' @return A list with components `exon` and `utr`.
#' @export
buildTranscriptAnnotation <- function(gtf) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Please install package 'rtracklayer'.")
  }

  if (!file.exists(gtf)) {
    stop("GTF file not found: ", gtf)
  }

  gr <- rtracklayer::import(gtf)
  df <- as.data.frame(gr)

  if (!"type" %in% colnames(df)) {
    stop("Imported GTF does not contain `type` column.")
  }

  std_df <- function(d) {
    out <- data.frame(
      chrom = as.character(d$seqnames),
      start = as.integer(d$start),
      end = as.integer(d$end),
      strand = as.character(d$strand),
      gene_id = if ("gene_id" %in% colnames(d)) as.character(d$gene_id) else NA_character_,
      gene_name = if ("gene_name" %in% colnames(d)) as.character(d$gene_name) else NA_character_,
      transcript_id = if ("transcript_id" %in% colnames(d)) as.character(d$transcript_id) else NA_character_,
      stringsAsFactors = FALSE
    )
    bad_gene_name <- is.na(out$gene_name) | out$gene_name == ""
    out$gene_name[bad_gene_name] <- out$gene_id[bad_gene_name]
    out
  }

  exon_df <- std_df(df[df$type == "exon", , drop = FALSE])

  utr_types <- c("UTR", "five_prime_UTR", "three_prime_UTR", "5UTR", "3UTR")
  utr_df <- std_df(df[df$type %in% utr_types, , drop = FALSE])

  list(
    exon = exon_df,
    utr = utr_df
  )
}
