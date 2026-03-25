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
.filter_chrX_by_sex <- function(df, sex = NULL,
                                chrX = c("X", "chrX"),
                                chrY = c("Y", "chrY"),
                                verbose = TRUE) {

  # ---- 先去掉 Y（默认行为）----
  df <- df[!df$contig %in% chrY, , drop = FALSE]

  # ---- 没传 sex → 去掉 X ----
  if (is.null(sex)) {
    if (isTRUE(verbose)) {
      warning("No sex information provided: chrX sites will be removed.")
    }
    return(df[!df$contig %in% chrX, , drop = FALSE])
  }

  # ---- 转换 sex 输入 ----
  if (is.data.frame(sex)) {
    if (!all(c("sample_id", "sex") %in% colnames(sex))) {
      stop("`sex` data.frame must contain columns: sample_id, sex")
    }
    sex_vec <- setNames(as.character(sex$sex), sex$sample_id)
  } else {
    sex_vec <- sex
  }

  raw_sex_vec <- sex_vec

  # ---- 标准化 ----
  sex_vec <- .normalize_sex(sex_vec)
  if (isTRUE(verbose)) print(table(sex_vec, useNA = "ifany"))

  # ---- 匹配到 df ----
  df$sex_tmp <- sex_vec[df$sample_id]

  # ---- 检查：未匹配样本 ----
  n_na_match <- sum(is.na(df$sex_tmp))
  if (n_na_match > 0 && isTRUE(verbose)) {
    warning(sprintf(
      "%d rows have no matching sex information (set to 'unknown').",
      n_na_match
    ))
  }

  df$sex_tmp[is.na(df$sex_tmp)] <- "unknown"

  # ---- 检查：标准化后分布 ----
  if (isTRUE(verbose)) {
    sex_table <- table(df$sex_tmp)

    if (length(sex_table) == 1 && names(sex_table) == "unknown") {
      warning("All samples are 'unknown' after sex normalization. chrX will be removed.")
    }

    if ("unknown" %in% names(sex_table)) {
      warning(sprintf(
        "%d entries have 'unknown' sex after normalization.",
        sex_table["unknown"]
      ))
    }

    if (!"female" %in% names(sex_table)) {
      warning("No female samples detected. All chrX sites will be removed.")
    }
  }

  # ---- X 过滤 ----
  df$is_chrX <- df$contig %in% chrX

  keep <- !(df$is_chrX & df$sex_tmp != "female")
  df <- df[keep, , drop = FALSE]

  # ---- 清理 ----
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
.normalize_sex <- function(sex_vec) {

  if (is.null(sex_vec)) return(NULL)

  nm <- names(sex_vec)

  sex_chr <- as.character(sex_vec)
  sex_chr <- trimws(tolower(sex_chr))

  female_set <- c("female", "f", "2", "xx")
  male_set   <- c("male", "m", "1", "xy")

  out <- rep("unknown", length(sex_chr))

  out[sex_chr %in% female_set] <- "female"
  out[sex_chr %in% male_set]   <- "male"

  names(out) <- nm
  out
}
