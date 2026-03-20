#' Run aberrant expression detection using OUTRIDER
#'
#' Perform cohort-based aberrant expression detection from either a featureCounts
#' table or a gene-by-sample raw count matrix.
#'
#' @param x A featureCounts table (data.frame with columns:
#' Geneid, Chr, Start, End, Strand, Length, followed by sample columns) or
#' a raw count matrix/data.frame.
#' @param geneLength Optional named numeric vector of gene lengths in bp.
#' Required when `x` is not a featureCounts table and `filterByFPKM = TRUE`.
#' @param colData Optional sample metadata as a data.frame.
#' @param rowData Optional gene annotation as a data.frame.
#' @param filterByFPKM Logical; whether to filter genes using FPKM.
#' @param minFPKM Numeric cutoff for FPKM filtering.
#' @param padjCutoff Adjusted p-value cutoff for defining outliers.
#' @param zCutoff Optional absolute z-score cutoff for defining outliers.
#' @param implementation Currently only "OUTRIDER".
#' @param returnModel Logical; whether to return the fitted OUTRIDER object.
#' @param ... Additional arguments passed to `OUTRIDER::OUTRIDER()`.
#'
#' @return An object of class `AEResult`.
#' @export
runAE <- function(
    x,
    geneLength = NULL,
    colData = NULL,
    rowData = NULL,
    filterByFPKM = TRUE,
    minFPKM = 1,
    padjCutoff = 0.05,
    zCutoff = NULL,
    implementation = "OUTRIDER",
    returnModel = FALSE,
    ...
) {
  .check_package_OUTRIDER()

  is_fc <- is.data.frame(x) &&
    all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% colnames(x))

  if (is_fc) {
    parsed <- parseFeatureCounts(x)
    counts <- parsed$counts
    geneLength <- parsed$geneLength
    if (is.null(rowData)) {
      rowData <- parsed$rowData
    }
  } else {
    counts <- .coerce_counts_matrix(x)
  }

  colData <- .coerce_coldata(colData, colnames(counts))
  rowData <- .coerce_rowdata(rowData, rownames(counts))

  if (!identical(implementation, "OUTRIDER")) {
    stop("Currently only implementation = 'OUTRIDER' is supported.")
  }

  fpkm <- NULL
  if (isTRUE(filterByFPKM)) {
    if (is.null(geneLength)) {
      stop("`geneLength` is required for FPKM filtering.")
    }
    fpkm <- computeFPKM(counts, geneLength)
  }

  keep <- rowSums(counts > 0, na.rm = TRUE) > 0
  if (isTRUE(filterByFPKM)) {
    keep <- keep & (rowSums(fpkm >= minFPKM, na.rm = TRUE) > 0)
  }

  counts_f <- counts[keep, , drop = FALSE]

  if (nrow(counts_f) == 0) {
    stop("No genes remaining after filtering.")
  }

  if (is.null(colData)) {
    colData_f <- data.frame(
      sampleID = colnames(counts_f),
      row.names = colnames(counts_f),
      stringsAsFactors = FALSE
    )
  } else {
    colData_f <- colData[colnames(counts_f), , drop = FALSE]
    if (!"sampleID" %in% colnames(colData_f)) {
      colData_f$sampleID <- rownames(colData_f)
    }
  }

  if (is.null(rowData)) {
    rowData_f <- data.frame(row.names = rownames(counts_f))
  } else {
    rowData_f <- rowData[rownames(counts_f), , drop = FALSE]
  }

  ods <- OUTRIDER::OutriderDataSet(
    countData = counts_f,
    colData = colData_f,
    rowData = rowData_f
  )

  ods <- OUTRIDER::filterExpression(
    ods,
    minCounts = TRUE,
    filterGenes = TRUE
  )

  ods <- OUTRIDER::OUTRIDER(ods, ...)
  res <- OUTRIDER::results(ods, all = TRUE)
  res <- as.data.frame(res, stringsAsFactors = FALSE)

  # 标准化列名
  res$sample_id <- res$sampleID
  res$gene_id <- res$geneID
  res$pvalue <- res$pValue
  res$padj <- res$padjust
  res$zscore <- res$zScore
  res$log2fc <- res$l2fc
  res$raw_count <- res$rawcounts
  res$norm_count <- res$normcounts
  res$expected_count <- res$normcounts

  # 加回 gene annotation
  if (!is.null(rowData_f) && nrow(res) > 0) {
    rd <- rowData_f
    rd$gene_id <- rownames(rd)
    idx <- match(res$gene_id, rd$gene_id)

    if ("gene_name" %in% colnames(rd)) {
      res$gene_name <- rd$gene_name[idx]
    }
    if ("fc_chr" %in% colnames(rd)) {
      res$chromosome <- rd$fc_chr[idx]
    }
  }

  # 定义 outlier
  res$is_outlier <- !is.na(res$padj) & res$padj < padjCutoff
  if (!is.null(zCutoff)) {
    res$is_outlier <- res$is_outlier &
      !is.na(res$zscore) &
      abs(res$zscore) >= zCutoff
  }

  # 保留标准列
  standard_cols <- c(
    "sample_id",
    "gene_id",
    "gene_name",
    "chromosome",
    "raw_count",
    "norm_count",
    "expected_count",
    "zscore",
    "log2fc",
    "pvalue",
    "padj",
    "is_outlier"
  )

  keep_cols <- intersect(standard_cols, colnames(res))
  res <- res[, keep_cols, drop = FALSE]

  model_out <- if (isTRUE(returnModel)) list(ods = ods) else NULL

  new_AEResult(
    result = res,
    metadata = list(
      implementation = implementation,
      n_genes_input = nrow(counts),
      n_genes_after_prefilter = nrow(counts_f),
      n_samples = ncol(counts),
      used_fpkm_filter = filterByFPKM,
      colData = colData_f
    ),
    parameters = list(
      minFPKM = minFPKM,
      padjCutoff = padjCutoff,
      zCutoff = zCutoff
    ),
    model = model_out
  )
}
