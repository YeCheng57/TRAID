#' Run aberrant expression detection using OUTRIDER
#'
#' Perform cohort-based aberrant expression detection for a target sample from
#' either a featureCounts table or a gene-by-sample raw count matrix.
#'
#' @param x A featureCounts table (data.frame with first 6 columns:
#' Geneid, Chr, Start, End, Strand, Length) or a raw count matrix/data.frame.
#' @param sampleID Character scalar specifying the target sample.
#' @param geneLength Optional named numeric vector of gene lengths in bp.
#' Required when `x` is not a featureCounts table and `filterByFPKM = TRUE`.
#' @param colData Optional sample metadata as a data.frame.
#' @param rowData Optional gene annotation as a data.frame.
#' @param filterByFPKM Logical; whether to filter genes using FPKM.
#' @param minFPKM Numeric cutoff for FPKM filtering.
#' @param padjCutoff Adjusted p-value cutoff.
#' @param zCutoff Optional absolute z-score cutoff.
#' @param implementation Currently only "OUTRIDER".
#' @param returnModel Logical; whether to return the fitted OUTRIDER object.
#' @param ... Additional arguments passed to `OUTRIDER::OUTRIDER()`.
#'
#' @return An object of class `AEResult`.
#' @export
runAE <- function(
    x,
    sampleID,
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

  .check_sample_id(sampleID, colnames(counts))
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
  res <- OUTRIDER::results(ods)
  res <- as.data.frame(res, stringsAsFactors = FALSE)

  if ("padjust" %in% colnames(res) && !"padj" %in% colnames(res)) {
    res$padj <- res$padjust
  }
  if ("pValue" %in% colnames(res) && !"pvalue" %in% colnames(res)) {
    res$pvalue <- res$pValue
  }
  if ("zScore" %in% colnames(res) && !"zscore" %in% colnames(res)) {
    res$zscore <- res$zScore
  }
  if ("sampleID" %in% colnames(res) && !"sample_id" %in% colnames(res)) {
    res$sample_id <- res$sampleID
  }
  if ("geneID" %in% colnames(res) && !"gene_id" %in% colnames(res)) {
    res$gene_id <- res$geneID
  }
  if ("rawcounts" %in% colnames(res) && !"raw_count" %in% colnames(res)) {
    res$raw_count <- res$rawcounts
  }
  if ("normcounts" %in% colnames(res) && !"expected_count" %in% colnames(res)) {
    res$expected_count <- res$normcounts
  }
  if ("l2fc" %in% colnames(res) && !"log2fc" %in% colnames(res)) {
    res$log2fc <- res$l2fc
  }

  sample_res <- res[res$sample_id == sampleID, , drop = FALSE]

  if (!is.null(rowData) && nrow(sample_res) > 0) {
    rd <- rowData_f
    rd$gene_id <- rownames(rd)
    idx <- match(sample_res$gene_id, rd$gene_id)
    extra_cols <- setdiff(colnames(rd), colnames(sample_res))
    for (nm in extra_cols) {
      sample_res[[nm]] <- rd[[nm]][idx]
    }
  }

  sample_res$is_outlier <- !is.na(sample_res$padj) & sample_res$padj < padjCutoff
  if (!is.null(zCutoff)) {
    sample_res$is_outlier <- sample_res$is_outlier &
      !is.na(sample_res$zscore) &
      abs(sample_res$zscore) >= zCutoff
  }

  preferred <- c(
    "sample_id", "gene_id", "gene_name",
    "raw_count", "expected_count",
    "zscore", "log2fc", "pvalue", "padj", "is_outlier"
  )
  keep_cols <- intersect(preferred, colnames(sample_res))
  other_cols <- setdiff(colnames(sample_res), keep_cols)
  sample_res <- sample_res[, c(keep_cols, other_cols), drop = FALSE]

  model_out <- if (isTRUE(returnModel)) list(ods = ods) else NULL

  new_AEResult(
    result = sample_res,
    metadata = list(
      sampleID = sampleID,
      implementation = implementation,
      n_genes_input = nrow(counts),
      n_genes_after_prefilter = nrow(counts_f),
      n_samples = ncol(counts),
      used_fpkm_filter = filterByFPKM
    ),
    parameters = list(
      minFPKM = minFPKM,
      padjCutoff = padjCutoff,
      zCutoff = zCutoff
    ),
    model = model_out
  )
}
