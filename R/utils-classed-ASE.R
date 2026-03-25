new_ASEResult <- function(result, metadata, parameters) {
  out <- list(
    result = result,
    metadata = metadata,
    parameters = parameters
  )
  class(out) <- c("ASEResult", "list")
  out
}

#' @export
print.ASEResult <- function(x, ...) {
  cat("<ASEResult>\n")
  cat(" Rows      :", nrow(x$result), "\n")
  cat(" Samples   :", x$metadata$n_samples, "\n")
  cat(" returnAll :", x$metadata$returnAll, "\n")
  invisible(x)
}
