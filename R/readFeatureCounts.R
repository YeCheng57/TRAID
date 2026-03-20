#' Parse featureCounts output table
#'
#' @param x A featureCounts result table as a data.frame.
#'
#' @return A list with counts matrix, gene length, and rowData.
#' @keywords internal
parseFeatureCounts <- function(x) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.")
  }

  req <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  miss <- setdiff(req, colnames(x))
  if (length(miss) > 0) {
    stop("Missing required featureCounts columns: ",
         paste(miss, collapse = ", "))
  }

  sample_cols <- setdiff(colnames(x), req)
  if (length(sample_cols) == 0) {
    stop("No sample count columns found in featureCounts table.")
  }

  counts <- as.matrix(x[, sample_cols, drop = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- x$Geneid

  geneLength <- x$Length
  names(geneLength) <- x$Geneid

  rowData <- data.frame(
    gene_id = x$Geneid,
    chr = x$Chr,
    start = x$Start,
    end = x$End,
    strand = x$Strand,
    length = x$Length,
    row.names = x$Geneid,
    stringsAsFactors = FALSE
  )

  list(
    counts = counts,
    geneLength = geneLength,
    rowData = rowData
  )
}
