.order_samples_for_plot <- function(
    df,
    sample_col = "sample_id",
    value_col,
    group_col = NULL,
    decreasing = TRUE
) {
  if (!sample_col %in% colnames(df)) {
    stop("`sample_col` not found in data.")
  }
  if (!value_col %in% colnames(df)) {
    stop("`value_col` not found in data.")
  }

  # 无分组：整体排序
  if (is.null(group_col) || !group_col %in% colnames(df) || all(is.na(df[[group_col]]))) {
    ord <- order(df[[value_col]], decreasing = decreasing, na.last = TRUE)
    df <- df[ord, , drop = FALSE]
    df[[sample_col]] <- factor(df[[sample_col]], levels = df[[sample_col]])
    return(df)
  }

  # 有分组：先组间，再组内排序
  group_chr <- as.character(df[[group_col]])

  if (decreasing) {
    ord <- order(group_chr, -xtfrm(df[[value_col]]), na.last = TRUE)
  } else {
    ord <- order(group_chr, xtfrm(df[[value_col]]), na.last = TRUE)
  }

  df <- df[ord, , drop = FALSE]
  df[[sample_col]] <- factor(df[[sample_col]], levels = df[[sample_col]])
  df
}
