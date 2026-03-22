#' Write exon annotation table
#'
#' @param x Exon annotation data.frame.
#' @param file Output file path.
#'
#' @return Invisibly returns `file`.
#' @export
writeExonAnnotation <- function(x, file) {
  .check_exon_annotation(x)

  utils::write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(file)
}
