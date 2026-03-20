.attach_group_info <- function(df, x, groupBy = NULL) {
  # 永远保证 group 列存在，且长度与 df 一致
  df$group <- rep(NA_character_, nrow(df))

  if (is.null(groupBy)) {
    return(df)
  }

  if (!inherits(x, "AEResult")) {
    return(df)
  }

  md <- x$metadata
  if (is.null(md) || is.null(md$colData)) {
    warning("`groupBy` was provided, but no `colData` is stored in the object. Plotting without grouping.")
    return(df)
  }

  colData <- md$colData
  if (!is.data.frame(colData)) {
    warning("Stored `colData` is not a data.frame. Plotting without grouping.")
    return(df)
  }

  if (!groupBy %in% colnames(colData)) {
    warning("`groupBy` column not found in `colData`. Plotting without grouping.")
    return(df)
  }

  if (is.null(rownames(colData))) {
    warning("`colData` has no row names. Plotting without grouping.")
    return(df)
  }

  idx <- match(df$sample_id, rownames(colData))
  df$group <- as.character(colData[[groupBy]][idx])
  df
}
