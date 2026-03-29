.attach_group_info <- function(df, x, groupBy = NULL, sample_info = NULL) {
  if (is.null(groupBy)) {
    return(df)
  }

  si <- NULL
  if (!is.null(sample_info)) {
    si <- sample_info
  } else if (!is.null(x$metadata$sample_info)) {
    si <- x$metadata$sample_info
  }

  if (is.null(si)) {
    warning("No sample_info available for grouping.")
    return(df)
  }

  if (!groupBy %in% colnames(si)) {
    warning("`groupBy` column not found in sample_info: ", groupBy)
    return(df)
  }

  if (!"sample_id" %in% colnames(df)) {
    return(df)
  }

  m <- match(df$sample_id, si$sample_id)
  df$group <- si[[groupBy]][m]

  if (any(is.na(df$group))) {
    warning("Some samples do not have group information.")
  }

  df
}
