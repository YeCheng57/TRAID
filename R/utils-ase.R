.annotateASERegion <- function(df, annotation) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE) ||
      !requireNamespace("IRanges", quietly = TRUE)) {
    stop("Please install packages 'GenomicRanges' and 'IRanges'.")
  }

  if (!is.list(annotation) || !all(c("exon", "utr") %in% names(annotation))) {
    stop("`annotation` must be a list with components `exon` and `utr`.")
  }

  # ---- 构建 query ----
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

  # ---- 初始化 ----
  df$region <- "other"
  df$gene_name <- NA_character_
  df$gene_id <- NA_character_

  # ---- exon overlaps ----
  if (length(exon_gr) > 0) {
    ov_exon <- GenomicRanges::findOverlaps(query, exon_gr)

    if (length(ov_exon) > 0) {
      q_hits <- S4Vectors::queryHits(ov_exon)
      s_hits <- S4Vectors::subjectHits(ov_exon)

      df$region[unique(q_hits)] <- "exonic"

      # gene annotation（允许多gene collapse）
      gene_map <- tapply(
        exon_df$gene_name[s_hits],
        q_hits,
        function(x) paste(unique(x[!is.na(x)]), collapse = ";")
      )

      df$gene_name[as.integer(names(gene_map))] <- gene_map

      gene_id_map <- tapply(
        exon_df$gene_id[s_hits],
        q_hits,
        function(x) paste(unique(x[!is.na(x)]), collapse = ";")
      )

      df$gene_id[as.integer(names(gene_id_map))] <- gene_id_map
    }
  }

  # ---- UTR overlaps（优先级高于 exon）----
  if (length(utr_gr) > 0) {
    ov_utr <- GenomicRanges::findOverlaps(query, utr_gr)

    if (length(ov_utr) > 0) {
      q_hits <- S4Vectors::queryHits(ov_utr)
      s_hits <- S4Vectors::subjectHits(ov_utr)

      df$region[unique(q_hits)] <- "UTR"

      gene_map <- tapply(
        utr_df$gene_name[s_hits],
        q_hits,
        function(x) paste(unique(x[!is.na(x)]), collapse = ";")
      )

      df$gene_name[as.integer(names(gene_map))] <- gene_map

      gene_id_map <- tapply(
        utr_df$gene_id[s_hits],
        q_hits,
        function(x) paste(unique(x[!is.na(x)]), collapse = ";")
      )

      df$gene_id[as.integer(names(gene_id_map))] <- gene_id_map
    }
  }

  df
}
.filter_chrX_by_sex <- function(df, sex = NULL, chrX = c("X", "chrX")) {

  # 没传 sex：默认去掉 X
  if (is.null(sex)) {
    return(df[!df$contig %in% chrX, , drop = FALSE])
  }

  # data.frame -> named vector
  if (is.data.frame(sex)) {
    if (!all(c("sample_id", "sex") %in% colnames(sex))) {
      stop("`sex` data.frame must contain columns: sample_id, sex")
    }
    sex_vec <- setNames(as.character(sex$sex), sex$sample_id)
  } else {
    sex_vec <- sex
  }

  sex_vec <- tolower(as.character(sex_vec))

  df$is_chrX <- df$contig %in% chrX
  df$sex_tmp <- sex_vec[df$sample_id]
  df$sex_tmp[is.na(df$sex_tmp)] <- "unknown"

  # 只有 female 保留 X
  keep <- !(df$is_chrX & df$sex_tmp != "female")
  df <- df[keep, , drop = FALSE]

  df$is_chrX <- NULL
  df$sex_tmp <- NULL

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

