new_AEResult <- function(result, metadata, parameters, model = NULL) {
  out <- list(
    result = result,
    metadata = metadata,
    parameters = parameters,
    model = model
  )
  class(out) <- c("AEResult", "list")
  out
}
#' @export
print.AEResult <- function(x, ...) {
  cat("<AEResult>\n")
  cat(" Sample:", x$metadata$sampleID, "\n")
  cat(" Genes :", nrow(x$result), "\n")
  if ("is_outlier" %in% colnames(x$result)) {
    cat(" Outliers:", sum(x$result$is_outlier, na.rm = TRUE), "\n")
  }
  invisible(x)
}
