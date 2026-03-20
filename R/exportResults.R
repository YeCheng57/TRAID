#' Export module results
#'
#' @param x A result object.
#' @param file Output file path.
#' @param which Component to export.
#'
#' @export
exportResults <- function(x, file, which = "result") {
  if (!which %in% names(x)) {
    stop("`which` not found in object.")
  }

  utils::write.table(
    x[[which]],
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(file)
}
