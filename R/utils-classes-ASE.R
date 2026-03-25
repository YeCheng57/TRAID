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

  if (!is.null(x$result) && nrow(x$result) > 0) {
    cat(" Rows            :", nrow(x$result), "\n")
    cat(" Samples         :", length(unique(x$result$sample_id)), "\n")
    cat(" Unique variants :", length(unique(x$result$variant_id)), "\n")
  } else {
    cat(" Rows            : 0\n")
  }

  cat(" ReturnAll       :", x$metadata$returnAll, "\n")

  cat("\nParameters:\n")
  cat(" minDP           :", x$parameters$minDP, "\n")
  cat(" pCutoff         :", x$parameters$pCutoff, "\n")
  cat(" ratioCutoff     :", x$parameters$ratioCutoff, "\n")

  if (!is.null(x$parameters$sex_provided)) {
    cat(" sex_provided    :", x$parameters$sex_provided, "\n")
  }

  invisible(x)
}
#' Summary of ASE results
#'
#' @param object An `ASEResult` object.
#' @param bySample Logical; whether to compute per-sample summary.
#'
#' @return A list containing summary statistics.
#' @export
summary.ASEResult <- function(object, bySample = TRUE, ...) {
  if (!inherits(object, "ASEResult")) {
    stop("`object` must be an ASEResult.")
  }

  df <- object$result

  if (nrow(df) == 0) {
    return(list(
      total_sites = 0,
      outlier_sites = 0
    ))
  }

  # ---- overall ----
  total_sites <- nrow(df)

  if ("is_outlier" %in% colnames(df)) {
    outlier_sites <- sum(df$is_outlier, na.rm = TRUE)
  } else {
    outlier_sites <- NA_integer_
  }

  region_table <- if ("region" %in% colnames(df)) {
    table(df$region)
  } else {
    NULL
  }

  gene_count <- if ("gene_name" %in% colnames(df)) {
    length(unique(df$gene_name[!is.na(df$gene_name)]))
  } else {
    NA_integer_
  }

  # ---- per sample ----
  sample_summary <- NULL

  if (isTRUE(bySample)) {
    sample_ids <- unique(df$sample_id)

    sample_summary <- data.frame(
      sample_id = sample_ids,
      n_sites = NA_integer_,
      n_outliers = NA_integer_,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(sample_ids)) {
      sid <- sample_ids[i]
      sub <- df[df$sample_id == sid, , drop = FALSE]

      sample_summary$n_sites[i] <- nrow(sub)

      if ("is_outlier" %in% colnames(sub)) {
        sample_summary$n_outliers[i] <- sum(sub$is_outlier, na.rm = TRUE)
      }
    }
  }

  out <- list(
    total_sites = total_sites,
    outlier_sites = outlier_sites,
    outlier_fraction = outlier_sites / total_sites,
    region_distribution = region_table,
    n_genes = gene_count,
    sample_summary = sample_summary
  )

  class(out) <- "summary.ASEResult"
  out
}
#' @export
print.summary.ASEResult <- function(x, ...) {

  cat("<ASE Summary>\n\n")

  cat("Total sites      :", x$total_sites, "\n")
  cat("Outlier sites    :", x$outlier_sites, "\n")

  if (!is.null(x$outlier_fraction) && !is.na(x$outlier_fraction)) {
    cat("Outlier fraction:", sprintf("%.4f", x$outlier_fraction), "\n")
  }

  if (!is.null(x$n_genes)) {
    cat("Genes involved   :", x$n_genes, "\n")
  }

  cat("\nRegion distribution:\n")
  if (!is.null(x$region_distribution)) {
    print(x$region_distribution)
  } else {
    cat("  (no region info)\n")
  }

  if (!is.null(x$sample_summary)) {
    cat("\nPer-sample summary:\n")
    print(x$sample_summary)
  }

  invisible(x)
}
