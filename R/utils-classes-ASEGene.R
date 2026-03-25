new_ASEGeneResult <- function(result, metadata, parameters) {
  out <- list(
    result = result,
    metadata = metadata,
    parameters = parameters
  )
  class(out) <- c("ASEGeneResult", "list")
  out
}

#' @export
print.ASEGeneResult <- function(x, ...) {
  cat("<ASEGeneResult>\n")

  if (!is.null(x$result) && nrow(x$result) > 0) {
    cat(" Rows         :", nrow(x$result), "\n")
    cat(" Samples      :", length(unique(x$result$sample_id)), "\n")
    cat(" Unique genes :", length(unique(x$result$gene_label)), "\n")

    if ("has_ASE" %in% colnames(x$result)) {
      cat(" ASE genes    :", sum(x$result$has_ASE, na.rm = TRUE), "\n")
    }
  } else {
    cat(" Rows         : 0\n")
  }

  cat(" ReturnAll    :", x$metadata$returnAll, "\n")

  cat("\nParameters:\n")
  cat(" useGene      :", x$parameters$useGene, "\n")
  cat(" minSites     :", x$parameters$minSitesPerGene, "\n")
  cat(" requireAll   :", x$parameters$requireAll, "\n")

  invisible(x)
}
