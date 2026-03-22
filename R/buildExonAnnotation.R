#' Build exon annotation from a local GTF file
#'
#' Parse a local GTF/GTF.GZ file and extract exon-level annotation into a
#' standardized data.frame.
#'
#' @param gtf Path to a local GTF or GTF.GZ file.
#'
#' @return A data.frame containing standardized exon annotation with columns:
#' \describe{
#'   \item{chrom}{Character. Chromosome name.}
#'   \item{start}{Integer. Exon start coordinate (1-based, inclusive).}
#'   \item{end}{Integer. Exon end coordinate (1-based, inclusive).}
#'   \item{strand}{Character. Strand information ("+", "-", or "*").}
#'   \item{gene_id}{Character. Gene identifier.}
#'   \item{gene_name}{Character. Gene symbol or gene name.}
#'   \item{transcript_id}{Character. Transcript identifier.}
#'   \item{exon_id}{Character. Unique exon identifier.}
#'   \item{exon_number}{Integer. Exon number within transcript, when available.}
#' }
#'
#' @details
#' This function only accepts a local file path. To support local files, URLs,
#' or cached loading in plotting functions, use the internal helper
#' \code{.resolve_exon_annotation()}.
#'
#' Chromosome naming must be consistent with the junction table used for
#' downstream plotting (e.g. both use `"1"` or both use `"chr1"`).
#'
#' @export
buildExonAnnotation <- function(gtf) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Please install package 'rtracklayer'.")
  }

  if (!is.character(gtf) || length(gtf) != 1L) {
    stop("`gtf` must be a character scalar.")
  }

  if (grepl("^https?://", gtf)) {
    stop("`buildExonAnnotation()` expects a local file path. For URL input, use a plotting function with `gtf=` or call `.resolve_exon_annotation()` internally.")
  }

  if (!file.exists(gtf)) {
    stop("GTF file not found: ", gtf)
  }

  message("Importing GTF with rtracklayer...")
  gr <- rtracklayer::import(gtf)

  gr <- gr[gr$type == "exon"]

  if (length(gr) == 0) {
    stop("No exon entries found in GTF.")
  }

  df <- as.data.frame(gr)

  exon_df <- data.frame(
    chrom = as.character(df$seqnames),
    start = as.integer(df$start),
    end = as.integer(df$end),
    strand = as.character(df$strand),
    gene_id = if ("gene_id" %in% colnames(df)) as.character(df$gene_id) else NA_character_,
    gene_name = if ("gene_name" %in% colnames(df)) as.character(df$gene_name) else NA_character_,
    transcript_id = if ("transcript_id" %in% colnames(df)) as.character(df$transcript_id) else NA_character_,
    exon_number = if ("exon_number" %in% colnames(df)) suppressWarnings(as.integer(df$exon_number)) else NA_integer_,
    stringsAsFactors = FALSE
  )

  bad_gene_name <- is.na(exon_df$gene_name) | exon_df$gene_name == ""
  exon_df$gene_name[bad_gene_name] <- exon_df$gene_id[bad_gene_name]

  bad_tx <- is.na(exon_df$transcript_id) | exon_df$transcript_id == ""
  exon_df$transcript_id[bad_tx] <- paste0("tx_", seq_len(sum(bad_tx)))

  exon_df$exon_id <- paste0(
    exon_df$gene_id, "_",
    exon_df$transcript_id, "_",
    exon_df$start, "_",
    exon_df$end
  )

  exon_df <- exon_df[, c(
    "chrom", "start", "end", "strand",
    "gene_id", "gene_name", "transcript_id",
    "exon_id", "exon_number"
  )]

  exon_df <- exon_df[order(
    exon_df$chrom,
    exon_df$gene_id,
    exon_df$transcript_id,
    exon_df$start,
    exon_df$end
  ), , drop = FALSE]

  rownames(exon_df) <- NULL

  .check_exon_annotation(exon_df)

  message("Done. Parsed ", nrow(exon_df), " exons.")
  exon_df
}
