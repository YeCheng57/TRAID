#' Find exon context for a focal junction
#'
#' @param chrom Chromosome.
#' @param start Junction start coordinate.
#' @param end Junction end coordinate.
#' @param gene_name Optional gene name used to restrict matching exons.
#' @param exon_annotation A standardized exon annotation data.frame.
#'
#' @return A list containing matched exon candidates and plotting coordinates.
#' @export
findJunctionExons <- function(
    chrom,
    start,
    end,
    gene_name = NULL,
    exon_annotation
) {
  .check_exon_annotation(exon_annotation)

  ex <- exon_annotation[exon_annotation$chrom == chrom, , drop = FALSE]

  if (!is.null(gene_name)) {
    ex2 <- ex[!is.na(ex$gene_name) & ex$gene_name == gene_name, , drop = FALSE]
    if (nrow(ex2) > 0) {
      ex <- ex2
    }
  }

  if (nrow(ex) == 0) {
    stop("No exon annotation found for requested chromosome/gene.")
  }

  # 命中 start/end 的 exon
  start_hits <- ex[ex$start <= start & ex$end >= start, , drop = FALSE]
  end_hits <- ex[ex$start <= end & ex$end >= end, , drop = FALSE]

  # 如果没有命中，找最近 exon
  nearest_exon <- function(pos, df) {
    mid <- (df$start + df$end) / 2
    df[which.min(abs(mid - pos)), , drop = FALSE]
  }

  if (nrow(start_hits) == 0) {
    start_hits <- nearest_exon(start, ex)
    start_hits$match_type <- "nearest"
  } else {
    start_hits$match_type <- "overlap"
  }

  if (nrow(end_hits) == 0) {
    end_hits <- nearest_exon(end, ex)
    end_hits$match_type <- "nearest"
  } else {
    end_hits$match_type <- "overlap"
  }

  # 选一个 transcript 优先级：
  # 1) start/end 同 transcript
  # 2) start exon 第一条
  chosen_tx <- NULL
  common_tx <- intersect(start_hits$transcript_id, end_hits$transcript_id)
  if (length(common_tx) > 0) {
    chosen_tx <- common_tx[1]
  } else if (nrow(start_hits) > 0) {
    chosen_tx <- start_hits$transcript_id[1]
  } else if (nrow(end_hits) > 0) {
    chosen_tx <- end_hits$transcript_id[1]
  }

  tx_exons <- ex[ex$transcript_id == chosen_tx, , drop = FALSE]
  tx_exons <- tx_exons[order(tx_exons$start, tx_exons$end), , drop = FALSE]

  start_exon <- start_hits[start_hits$transcript_id == chosen_tx, , drop = FALSE]
  if (nrow(start_exon) == 0) {
    start_exon <- start_hits[1, , drop = FALSE]
  }

  end_exon <- end_hits[end_hits$transcript_id == chosen_tx, , drop = FALSE]
  if (nrow(end_exon) == 0) {
    end_exon <- end_hits[1, , drop = FALSE]
  }

  list(
    transcript_id = chosen_tx,
    transcript_exons = tx_exons,
    start_exon = start_exon[1, , drop = FALSE],
    end_exon = end_exon[1, , drop = FALSE],
    start_hits = start_hits,
    end_hits = end_hits
  )
}
