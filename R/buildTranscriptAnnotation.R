#' Parse GATK ASEReadCounter table
#'
#' @param x A data.frame containing GATK ASEReadCounter results.
#' @param sampleID Optional sample ID. If `x` does not contain a sample column,
#' this value will be used.
#'
#' @return A standardized ASE table.
#' @export
parseASEReadCounter <- function(x, sampleID = NULL) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.")
  }

  cn <- colnames(x)

  # 常见列名兼容
  contig_col <- intersect(c("contig", "Contig", "CHROM", "chr"), cn)
  pos_col <- intersect(c("position", "Position", "POS", "pos"), cn)
  ref_col <- intersect(c("refAllele", "RefAllele", "REF", "ref"), cn)
  alt_col <- intersect(c("altAllele", "AltAllele", "ALT", "alt"), cn)
  ref_count_col <- intersect(c("refCount", "RefCount", "ref_count"), cn)
  alt_count_col <- intersect(c("altCount", "AltCount", "alt_count"), cn)
  total_col <- intersect(c("totalCount", "TotalCount", "DP", "dp", "total_count"), cn)
  sample_col <- intersect(c("sample_id", "sampleID", "Sample", "sample"), cn)

  reqs <- c(contig_col[1], pos_col[1], ref_col[1], alt_col[1], ref_count_col[1], alt_count_col[1])
  if (any(is.na(reqs)) || any(reqs == "")) {
    stop("Could not identify required ASEReadCounter columns.")
  }

  out <- data.frame(
    sample_id = if (length(sample_col) > 0) as.character(x[[sample_col[1]]]) else sampleID,
    contig = as.character(x[[contig_col[1]]]),
    position = as.integer(x[[pos_col[1]]]),
    ref = as.character(x[[ref_col[1]]]),
    alt = as.character(x[[alt_col[1]]]),
    ref_count = as.numeric(x[[ref_count_col[1]]]),
    alt_count = as.numeric(x[[alt_count_col[1]]]),
    stringsAsFactors = FALSE
  )

  if (is.null(out$sample_id)) {
    stop("No sample column found and `sampleID` was not provided.")
  }

  if (length(total_col) > 0) {
    out$total_count <- as.numeric(x[[total_col[1]]])
  } else {
    out$total_count <- out$ref_count + out$alt_count
  }

  # variant id
  out$variant_id <- paste(out$contig, out$position, out$ref, out$alt, sep = ":")

  # 默认先标成 SNV / non-SNV
  out$variant_type <- ifelse(
    nchar(out$ref) == 1 & nchar(out$alt) == 1,
    "SNV",
    "non-SNV"
  )

  # 如果原表已有 gene/region 信息，也带上
  if ("gene_name" %in% cn) out$gene_name <- as.character(x$gene_name)
  if ("gene_id" %in% cn) out$gene_id <- as.character(x$gene_id)
  if ("region" %in% cn) out$region <- as.character(x$region)

  out
}
