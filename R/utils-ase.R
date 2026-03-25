.annotateASERegion <- function(df, annotation) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop("Please install packages 'GenomicRanges' and 'IRanges'.")
  }

  if (!is.list(annotation) || !all(c("exon", "utr") %in% names(annotation))) {
    stop("`annotation` must be a list with components `exon` and `utr`.")
  }

  query <- GenomicRanges::GRanges(
    seqnames = df$contig,
    ranges = IRanges::IRanges(start = df$position, end = df$position)
  )

  exon_df <- annotation$exon
  utr_df <- annotation$utr

  exon_gr <- GenomicRanges::GRanges(
    seqnames = exon_df$chrom,
    ranges = IRanges::IRanges(start = exon_df$start, end = exon_df$end)
  )
  utr_gr <- GenomicRanges::GRanges(
    seqnames = utr_df$chrom,
    ranges = IRanges::IRanges(start = utr_df$start, end = utr_df$end)
  )

  df$region <- "other"

  if (length(exon_gr) > 0) {
    ov_exon <- GenomicRanges::findOverlaps(query, exon_gr)
    if (length(ov_exon) > 0) {
      df$region[unique(S4Vectors::queryHits(ov_exon))] <- "exonic"
    }
  }

  if (length(utr_gr) > 0) {
    ov_utr <- GenomicRanges::findOverlaps(query, utr_gr)
    if (length(ov_utr) > 0) {
      df$region[unique(S4Vectors::queryHits(ov_utr))] <- "UTR"
    }
  }

  df
}


.run_binom_ase <- function(df, pAdjustMethod = "fdr") {
  if (nrow(df) == 0) {
    df$pvalue <- numeric(0)
    df$padj <- numeric(0)
    df$alt_ratio <- numeric(0)
    df$major_ratio <- numeric(0)
    return(df)
  }

  df$alt_ratio <- df$alt_count / df$total_count
  df$major_ratio <- pmax(df$ref_count, df$alt_count) / df$total_count

  df$pvalue <- vapply(
    seq_len(nrow(df)),
    function(i) {
      stats::binom.test(
        x = df$alt_count[i],
        n = df$total_count[i],
        p = 0.5
      )$p.value
    },
    numeric(1)
  )

  df$padj <- stats::p.adjust(df$pvalue, method = pAdjustMethod)
  df
}
