.coerce_counts_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)) {
    stop("`counts` must be a matrix or data.frame.")
  }

  storage.mode(x) <- "numeric"

  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("`counts` must have row names (genes) and column names (samples).")
  }

  if (anyDuplicated(rownames(x))) {
    stop("Duplicate gene IDs found in row names of `counts`.")
  }
  if (anyDuplicated(colnames(x))) {
    stop("Duplicate sample IDs found in column names of `counts`.")
  }

  x
}

.coerce_coldata <- function(colData, sample_ids) {
  if (is.null(colData)) {
    return(NULL)
  }
  if (!is.data.frame(colData)) {
    stop("`colData` must be a data.frame.")
  }

  if (is.null(rownames(colData))) {
    if ("sample_id" %in% colnames(colData)) {
      rownames(colData) <- colData$sample_id
    } else if ("sampleID" %in% colnames(colData)) {
      rownames(colData) <- colData$sampleID
    } else {
      stop("`colData` must have row names or a sample_id/sampleID column.")
    }
  }

  miss <- setdiff(sample_ids, rownames(colData))
  if (length(miss) > 0) {
    stop("Missing samples in `colData`: ", paste(miss, collapse = ", "))
  }

  colData[sample_ids, , drop = FALSE]
}

.coerce_rowdata <- function(rowData, gene_ids) {
  if (is.null(rowData)) {
    return(NULL)
  }
  if (!is.data.frame(rowData)) {
    stop("`rowData` must be a data.frame.")
  }

  if (is.null(rownames(rowData))) {
    if ("gene_id" %in% colnames(rowData)) {
      rownames(rowData) <- rowData$gene_id
    } else if ("geneID" %in% colnames(rowData)) {
      rownames(rowData) <- rowData$geneID
    } else {
      stop("`rowData` must have row names or a gene_id/geneID column.")
    }
  }

  miss <- setdiff(gene_ids, rownames(rowData))
  if (length(miss) > 0) {
    stop("Missing genes in `rowData`: first few missing = ",
         paste(utils::head(miss, 10), collapse = ", "))
  }

  rowData[gene_ids, , drop = FALSE]
}

.check_package_OUTRIDER <- function() {
  if (!requireNamespace("OUTRIDER", quietly = TRUE)) {
    stop(
      "Package 'OUTRIDER' is required but not installed. ",
      "Please install it with BiocManager::install('OUTRIDER')."
    )
  }
  invisible(TRUE)
}
