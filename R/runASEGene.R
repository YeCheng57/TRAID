#' Run gene-level ASE analysis
#'
#' A preliminary implementation of gene-level ASE aggregation based on site-level
#' ASE results. A gene is considered to exhibit ASE when all retained
#' heterozygous sites within the gene satisfy the site-level ASE criterion.
#'
#' @param x An `ASEResult` object or a data.frame containing site-level ASE results.
#' @param useGene One of `"gene_name"` or `"gene_id"`.
#' @param minSitesPerGene Minimum number of heterozygous sites required per gene.
#' @param requireAll Logical; if TRUE, all retained sites in a gene must show ASE.
#' If FALSE, at least one retained site showing ASE is sufficient.
#' @param returnAll Logical; whether to return all tested genes. If FALSE, only genes with ASE are returned.
#'
#' @return An object of class `ASEGeneResult`.
#' @export
runASEGene <- function(
    x,
    useGene = c("gene_name", "gene_id"),
    minSitesPerGene = 2,
    requireAll = TRUE,
    returnAll = FALSE
) {
  useGene <- match.arg(useGene)

  if (inherits(x, "ASEResult")) {
    df <- x$result
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("`x` must be an ASEResult object or a data.frame.")
  }

  if (nrow(df) == 0) {
    return(new_ASEGeneResult(
      result = data.frame(),
      metadata = list(returnAll = returnAll),
      parameters = list(
        useGene = useGene,
        minSitesPerGene = minSitesPerGene,
        requireAll = requireAll
      )
    ))
  }

  if (!useGene %in% colnames(df)) {
    stop("Column `", useGene, "` not found in input.")
  }

  if (!"sample_id" %in% colnames(df)) {
    stop("Column `sample_id` not found in input.")
  }

  if (!"is_outlier" %in% colnames(df)) {
    stop("Column `is_outlier` not found in input.")
  }

  # 去掉没有 gene 标注的位点
  gene_label <- as.character(df[[useGene]])
  keep <- !is.na(gene_label) & gene_label != "" & gene_label != "NA"
  df <- df[keep, , drop = FALSE]

  if (nrow(df) == 0) {
    return(new_ASEGeneResult(
      result = data.frame(),
      metadata = list(returnAll = returnAll),
      parameters = list(
        useGene = useGene,
        minSitesPerGene = minSitesPerGene,
        requireAll = requireAll
      )
    ))
  }

  df$gene_label <- as.character(df[[useGene]])

  split_key <- paste(df$sample_id, df$gene_label, sep = "||")
  idx_list <- split(seq_len(nrow(df)), split_key)

  res_list <- lapply(idx_list, function(idx) {
    d <- df[idx, , drop = FALSE]

    n_sites <- nrow(d)
    n_ASE_sites <- sum(d$is_outlier, na.rm = TRUE)

    if (isTRUE(requireAll)) {
      gene_has_ASE <- n_sites >= minSitesPerGene && all(d$is_outlier)
    } else {
      gene_has_ASE <- n_sites >= minSitesPerGene && any(d$is_outlier)
    }

    data.frame(
      sample_id = d$sample_id[1],
      gene_label = d$gene_label[1],
      n_sites = n_sites,
      n_ASE_sites = n_ASE_sites,
      ASE_fraction = n_ASE_sites / n_sites,
      mean_major_ratio = if ("major_ratio" %in% colnames(d)) mean(d$major_ratio, na.rm = TRUE) else NA_real_,
      mean_alt_ratio = if ("alt_ratio" %in% colnames(d)) mean(d$alt_ratio, na.rm = TRUE) else NA_real_,
      has_ASE = gene_has_ASE,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res_list)
  rownames(out) <- NULL

  if (!isTRUE(returnAll)) {
    out <- out[!is.na(out$has_ASE) & out$has_ASE, , drop = FALSE]
    rownames(out) <- NULL
  }

  new_ASEGeneResult(
    result = out,
    metadata = list(
      returnAll = returnAll
    ),
    parameters = list(
      useGene = useGene,
      minSitesPerGene = minSitesPerGene,
      requireAll = requireAll
    )
  )
}
