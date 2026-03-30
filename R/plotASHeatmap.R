#' Plot AS heatmap across samples
#'
#' Visualize global splicing patterns across samples using PSI or delta PSI.
#'
#' @param x An `ASResult` object.
#' @param site One of `"start"`, `"end"`, or `"both"`.
#' @param value One of `"delta_psi"` or `"psi"`.
#' @param topN Number of top variable events to display.
#' @param sample_info Optional data.frame; overrides `x$metadata$sample_info`.
#' @param annotation_cols Optional character vector of sample_info columns.
#' @param cluster_rows Logical.
#' @param cluster_cols Logical.
#' @param scale One of `"row"` or `"none"`.
#' @param show_rownames Logical.
#' @param show_colnames Logical.
#'
#' @return Invisibly returns matrix and annotation.
#' @export
plotASHeatmap <- function(
    x,
    site = c("start", "end", "both"),
    value = c("delta_psi", "psi"),
    topN = 100,
    sample_info = NULL,
    annotation_cols = NULL,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    show_rownames = FALSE,
    show_colnames = FALSE
) {

  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Please install `reshape2`.")
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install `pheatmap`.")
  }

  site <- match.arg(site)
  value <- match.arg(value)

  # ---- 取数据 ----
  df <- switch(
    site,
    start = x$result_start,
    end = x$result_end,
    both = rbind(x$result_start, x$result_end)
  )

  if (is.null(df) || nrow(df) == 0) {
    stop("No AS results available.")
  }

  req <- c("sample_id", "chrom", "start", "end", value)
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  # ---- event_id（稳定唯一）----
  df$event_id <- paste0(
    ifelse(site == "both",
           paste0(ifelse(df$start_id == df$end_id, "event", "event")),
           site),
    "|",
    df$chrom, ":", df$start, "-", df$end
  )

  # ---- 保留必要列 ----
  df <- df[, c("event_id", "sample_id", value), drop = FALSE]

  # ---- 转矩阵 ----
  mat_df <- reshape2::dcast(df, event_id ~ sample_id, value.var = value)
  rownames(mat_df) <- mat_df$event_id
  mat_df$event_id <- NULL
  mat <- as.matrix(mat_df)

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Heatmap matrix is empty.")
  }

  storage.mode(mat) <- "numeric"

  # =========================================================
  # 🔥 核心修复（AS逻辑）
  # =========================================================

  # NA = 无junction → 0
  mat[is.na(mat)] <- 0

  # 去掉非有限值
  keep <- apply(mat, 1, function(z) all(is.finite(z)))
  mat <- mat[keep, , drop = FALSE]

  if (nrow(mat) == 0) {
    stop("No valid rows after filtering.")
  }

  # 去掉零方差（防止scale炸）
  rv <- apply(mat, 1, stats::var)
  keep <- is.finite(rv) & rv > 0
  mat <- mat[keep, , drop = FALSE]
  rv <- rv[keep]

  if (nrow(mat) == 0) {
    stop("No variable events available.")
  }

  # 选 topN
  keep_n <- min(topN, nrow(mat))
  idx <- order(rv, decreasing = TRUE)[seq_len(keep_n)]
  mat <- mat[idx, , drop = FALSE]

  # =========================================================
  # sample annotation
  # =========================================================

  si <- if (!is.null(sample_info)) sample_info else x$metadata$sample_info
  ann_col <- NULL

  if (!is.null(annotation_cols)) {

    if (is.null(si)) {
      warning("No sample_info available.")
    } else {

      if (!"sample_id" %in% colnames(si)) {
        stop("sample_info must contain `sample_id`.")
      }

      bad <- setdiff(annotation_cols, colnames(si))
      if (length(bad) > 0) {
        stop("Invalid annotation columns: ", paste(bad, collapse = ", "))
      }

      m <- match(colnames(mat), si$sample_id)

      ann_col <- si[m, annotation_cols, drop = FALSE]
      rownames(ann_col) <- si$sample_id[m]
    }
  }

  # =========================================================
  # 画图
  # =========================================================

  pheatmap::pheatmap(
    mat,
    annotation_col = ann_col,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale = scale,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    main = paste0("AS heatmap (", value, ", ", site, ")"),
    na_col = "grey90"
  )

  invisible(list(
    matrix = mat,
    annotation_col = ann_col
  ))
}
