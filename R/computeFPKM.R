#' Compute FPKM from count matrix
#'
#' @param counts Gene by sample raw count matrix.
#' @param geneLength Named numeric vector of gene lengths in bp.
#'
#' @return A numeric FPKM matrix.
#' @export
computeFPKM <- function(counts, geneLength) {
  counts <- as.matrix(counts)

  if (is.null(rownames(counts))) {
    stop("`counts` must have gene IDs as row names.")
  }
  if (is.null(names(geneLength))) {
    stop("`geneLength` must be a named vector with gene IDs.")
  }

  geneLength <- geneLength[rownames(counts)]

  if (any(is.na(geneLength))) {
    stop("Missing gene length for some genes.")
  }

  libSize <- colSums(counts, na.rm = TRUE)

  fpkm <- sweep(counts, 2, libSize, "/")
  fpkm <- sweep(fpkm, 1, geneLength, "/")
  fpkm <- fpkm * 1e9

  fpkm
}
