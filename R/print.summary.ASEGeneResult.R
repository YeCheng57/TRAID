#' @export
print.summary.ASEGeneResult <- function(x, ...) {

  cat("<ASE Gene-level Summary>\n\n")

  cat("Total genes         :", x$total_genes, "\n")
  cat("Genes with ASE      :", x$ASE_genes, "\n")

  if (!is.null(x$ASE_fraction) && !is.na(x$ASE_fraction)) {
    cat("ASE gene fraction   :", sprintf("%.4f", x$ASE_fraction), "\n")
  }

  if (!is.null(x$mean_sites_per_gene)) {
    cat("Mean sites per gene :", sprintf("%.2f", x$mean_sites_per_gene), "\n")
  }

  if (!is.null(x$mean_ASE_fraction_per_gene)) {
    cat("Mean ASE fraction   :", sprintf("%.4f", x$mean_ASE_fraction_per_gene), "\n")
  }

  if (!is.null(x$sample_summary)) {
    cat("\nPer-sample summary:\n")
    print(x$sample_summary)
  }

  invisible(x)
}
