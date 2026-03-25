#' Summary of gene-level ASE results
#'
#' @param object An `ASEGeneResult` object.
#' @param bySample Logical; whether to compute per-sample summary.
#'
#' @return A list containing summary statistics.
#' @export
summary.ASEGeneResult <- function(object, bySample = TRUE, ...) {
  if (!inherits(object, "ASEGeneResult")) {
    stop("`object` must be an ASEGeneResult.")
  }

  df <- object$result

  if (nrow(df) == 0) {
    return(list(
      total_genes = 0,
      ASE_genes = 0
    ))
  }

  # ---- overall ----
  total_genes <- nrow(df)

  ASE_genes <- if ("has_ASE" %in% colnames(df)) {
    sum(df$has_ASE, na.rm = TRUE)
  } else {
    NA_integer_
  }

  ASE_fraction <- ASE_genes / total_genes

  # ---- sites info ----
  mean_sites <- if ("n_sites" %in% colnames(df)) {
    mean(df$n_sites, na.rm = TRUE)
  } else {
    NA_real_
  }

  mean_ASE_fraction <- if ("ASE_fraction" %in% colnames(df)) {
    mean(df$ASE_fraction, na.rm = TRUE)
  } else {
    NA_real_
  }

  # ---- per sample ----
  sample_summary <- NULL

  if (isTRUE(bySample)) {
    sample_ids <- unique(df$sample_id)

    sample_summary <- data.frame(
      sample_id = sample_ids,
      n_genes = NA_integer_,
      n_ASE_genes = NA_integer_,
      ASE_fraction = NA_real_,
      stringsAsFactors = FALSE
    )

    for (i in seq_along(sample_ids)) {
      sid <- sample_ids[i]
      sub <- df[df$sample_id == sid, , drop = FALSE]

      sample_summary$n_genes[i] <- nrow(sub)

      if ("has_ASE" %in% colnames(sub)) {
        sample_summary$n_ASE_genes[i] <- sum(sub$has_ASE, na.rm = TRUE)
        sample_summary$ASE_fraction[i] <- sample_summary$n_ASE_genes[i] / sample_summary$n_genes[i]
      }
    }
  }

  out <- list(
    total_genes = total_genes,
    ASE_genes = ASE_genes,
    ASE_fraction = ASE_fraction,
    mean_sites_per_gene = mean_sites,
    mean_ASE_fraction_per_gene = mean_ASE_fraction,
    sample_summary = sample_summary
  )

  class(out) <- "summary.ASEGeneResult"
  out
}
