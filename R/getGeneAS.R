#' Extract AS results for one gene
#'
#' @param x An `ASResult` object.
#' @param geneID Gene identifier or gene symbol.
#' @param site One of `"start"`, `"end"`, or `"both"`.
#' @param use Column used for matching, usually `"gene_names"` or `"gene_ids"`.
#'
#' @return A data.frame if one site is requested, or a list if `site = "both"`.
#' @export
getGeneAS <- function(
    x,
    geneID,
    site = c("both","start", "end"),
    use = "gene_names"
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  use <- match.arg(use)

  extract_one <- function(df, geneID, use) {
    if (is.null(df)) {
      return(NULL)
    }
    if (!use %in% colnames(df)) {
      stop("Column `", use, "` not found in AS result.")
    }
    df[df[[use]] == geneID&df$is_outlier, , drop = FALSE]
  }

  if (site == "start") {
    return(extract_one(x$result_start, geneID, use))
  }
  if (site == "end") {
    return(extract_one(x$result_end, geneID, use))
  }

  rbind(
      start = extract_one(x$result_start, geneID, use),
      end = extract_one(x$result_end, geneID, use)
    )
}
