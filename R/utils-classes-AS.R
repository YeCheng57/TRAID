new_ASResult <- function(result_start, result_end, metadata, parameters) {
  out <- list(
    result_start = result_start,
    result_end = result_end,
    metadata = metadata,
    parameters = parameters
  )
  class(out) <- c("ASResult", "list")
  out
}

#' @export
print.ASResult <- function(x, ...) {
  cat("<ASResult>\n")
  cat(" Mode      :", x$metadata$mode, "\n")
  cat(" returnAll :", x$metadata$returnAll, "\n")

  if (!is.null(x$result_start)) {
    cat(" Start rows :", nrow(x$result_start), "\n")
  }
  if (!is.null(x$result_end)) {
    cat(" End rows   :", nrow(x$result_end), "\n")
  }

  invisible(x)
}
