#' Read exon annotation table
#'
#' @param file Path to a tab-delimited exon annotation file.
#'
#' @return A standardized exon annotation data.frame.
#' @export
readExonAnnotation <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  df <- utils::read.delim(
    file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  .check_exon_annotation(df)
  df
}
