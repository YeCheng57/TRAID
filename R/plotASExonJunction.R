#' Plot exon-aware junction cartoon for one AS event
#'
#' @param x An `ASResult` object, preferably generated with `returnAll = TRUE`.
#' @param sampleID Character scalar specifying the target sample.
#' @param chrom Chromosome of the focal junction.
#' @param start Junction start coordinate.
#' @param end Junction end coordinate.
#' @param site One of "start" or "end".
#' @param exon_annotation A standardized exon annotation data.frame.
#' @param gtf Optional local path or URL to a GTF/GTF.GZ file. Used only if
#'   `exon_annotation` is not provided.
#' @param cache Logical; whether to cache parsed exon annotation when `gtf` is used.
#' @param topN Integer; number of competing junctions to display in addition to the focal junction.
#' @param showReference Logical; whether to draw a reference transcript track.
#'
#' @return A ggplot object.
#' @export
plotASExonJunction <- function(
    x,
    sampleID,
    chrom,
    start,
    end,
    site = c("start", "end"),
    exon_annotation = NULL,
    gtf = NULL,
    cache = TRUE,
    topN = 2,
    showReference = TRUE
) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)
  df <- if (site == "start") x$result_start else x$result_end

  if (is.null(df) || nrow(df) == 0) {
    stop("Requested site result is empty.")
  }
  if (!sampleID %in% df$sample_id) {
    stop("`sampleID` not found in result.")
  }

  exons <- .resolve_exon_annotation(
    exon_annotation = exon_annotation,
    gtf = gtf,
    cache = cache
  )

  if (site == "start") {
    event_df <- df[df$chrom == chrom & df$start == start, , drop = FALSE]
    partner_col <- "end"
    focal_partner <- end
    anchor_pos <- start
  } else {
    event_df <- df[df$chrom == chrom & df$end == end, , drop = FALSE]
    partner_col <- "start"
    focal_partner <- start
    anchor_pos <- end
  }

  if (nrow(event_df) == 0) {
    stop("Requested focal event not found in AS result.")
  }

  gene_name <- NULL
  if ("gene_names" %in% colnames(event_df)) {
    gn <- unique(event_df$gene_names)
    gn <- gn[!is.na(gn) & gn != "" & gn != "NA"]
    if (length(gn) > 0) {
      gene_name <- gn[1]
    }
  }

  exon_ctx <- findJunctionExons(
    chrom = chrom,
    start = start,
    end = end,
    gene_name = gene_name,
    exon_annotation = exons
  )

  tx_exons_all <- exon_ctx$transcript_exons
  if (nrow(tx_exons_all) == 0) {
    stop("No transcript exons found for plotting.")
  }

  sample_df <- event_df[event_df$sample_id == sampleID, , drop = FALSE]
  others_df <- event_df[event_df$sample_id != sampleID, , drop = FALSE]

  if (nrow(sample_df) == 0) {
    stop("Target sample has no junctions at this focal site.")
  }

  # focal junction + 最主要竞争junction
  sample_rank <- sample_df[order(sample_df$junction_count, decreasing = TRUE), , drop = FALSE]
  keep_partners <- unique(c(focal_partner, head(sample_rank[[partner_col]], topN + 1L)))
  keep_partners <- keep_partners[!is.na(keep_partners)]

  sample_df <- sample_df[sample_df[[partner_col]] %in% keep_partners, , drop = FALSE]
  others_df <- others_df[others_df[[partner_col]] %in% keep_partners, , drop = FALSE]

  locate_exon_for_pos <- function(pos, ex_df) {
    hit <- ex_df[ex_df$start <= pos & ex_df$end >= pos, , drop = FALSE]
    if (nrow(hit) > 0) {
      return(hit[1, , drop = FALSE])
    }
    mid <- (ex_df$start + ex_df$end) / 2
    ex_df[which.min(abs(mid - pos)), , drop = FALSE]
  }

  match_exon_idx <- function(exon_row, tx) {
    idx <- which(tx$start == exon_row$start & tx$end == exon_row$end)
    if (length(idx) == 0) {
      which.min(abs(((tx$start + tx$end) / 2) - ((exon_row$start + exon_row$end) / 2)))
    } else {
      idx[1]
    }
  }

  anchor_exon <- locate_exon_for_pos(anchor_pos, tx_exons_all)

  partner_exons <- lapply(keep_partners, function(pos) {
    ex <- locate_exon_for_pos(pos, tx_exons_all)
    ex$partner_pos <- pos
    ex
  })
  partner_exons <- do.call(rbind, partner_exons)
  rownames(partner_exons) <- NULL
  partner_exons$is_focal <- partner_exons$partner_pos == focal_partner

  # 避免 merge 后 start/end 冲突
  partner_exons$exon_start <- partner_exons$start
  partner_exons$exon_end <- partner_exons$end

  # 局部 transcript 区域
  tx_ref_all <- tx_exons_all[order(tx_exons_all$start, tx_exons_all$end), , drop = FALSE]
  idx_anchor <- match_exon_idx(anchor_exon, tx_ref_all)
  idx_partners <- vapply(
    seq_len(nrow(partner_exons)),
    function(i) match_exon_idx(partner_exons[i, , drop = FALSE], tx_ref_all),
    numeric(1)
  )

  idx_min <- max(1, min(c(idx_anchor, idx_partners)) - 1)
  idx_max <- min(nrow(tx_ref_all), max(c(idx_anchor, idx_partners)) + 1)

  tx_ref <- tx_ref_all[idx_min:idx_max, , drop = FALSE]
  tx_ref <- tx_ref[order(tx_ref$start, tx_ref$end), , drop = FALSE]

  # 统一布局
  x_all <- seq_len(nrow(tx_ref)) * 2.4 + 2.0
  tx_ref$xmid <- x_all

  anchor_local_idx <- match_exon_idx(anchor_exon, tx_ref)
  x_anchor <- tx_ref$xmid[anchor_local_idx]

  partner_exons$xmid <- vapply(
    seq_len(nrow(partner_exons)),
    function(i) {
      idx <- match_exon_idx(partner_exons[i, , drop = FALSE], tx_ref)
      tx_ref$xmid[idx]
    },
    numeric(1)
  )

  exon_w <- 0.55
  exon_h <- 0.20
  y_sample <- 2.0
  y_ctrl <- 0.0
  y_ref <- -1.2

  map_pos_to_exon_x <- function(pos, xm, exon_start, exon_end, exon_w) {
    vals <- c(pos, xm, exon_start, exon_end)

    if (anyNA(vals)) {
      return(as.numeric(xm))
    }

    delta <- exon_end - exon_start
    if (!is.finite(delta) || delta <= 0) {
      return(as.numeric(xm))
    }

    rel <- (pos - exon_start) / delta
    rel <- max(0, min(1, rel))

    as.numeric(xm) - exon_w + rel * (2 * exon_w)
  }

  # controls mean ± sd
  if (nrow(others_df) > 0) {
    partner_vec <- others_df[[partner_col]]

    mean_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = partner_vec),
      FUN = mean,
      na.rm = TRUE
    )
    sd_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = partner_vec),
      FUN = stats::sd,
      na.rm = TRUE
    )

    colnames(mean_df)[2] <- "mean"
    colnames(sd_df)[2] <- "sd"

    others_sum <- merge(mean_df, sd_df, by = "partner", all = TRUE, sort = FALSE)
    others_sum$sd[is.na(others_sum$sd)] <- 0
  } else {
    others_sum <- data.frame(
      partner = numeric(0),
      mean = numeric(0),
      sd = numeric(0)
    )
  }

  # exon rectangles
  rect_df <- rbind(
    data.frame(
      panel = c("Sample", "Controls"),
      xmin = x_anchor - exon_w,
      xmax = x_anchor + exon_w,
      ymin = c(y_sample - exon_h / 2, y_ctrl - exon_h / 2),
      ymax = c(y_sample + exon_h / 2, y_ctrl + exon_h / 2),
      exon_type = "anchor",
      stringsAsFactors = FALSE
    ),
    data.frame(
      panel = rep(c("Sample", "Controls"), each = nrow(partner_exons)),
      xmin = rep(partner_exons$xmid - exon_w, 2),
      xmax = rep(partner_exons$xmid + exon_w, 2),
      ymin = c(rep(y_sample - exon_h / 2, nrow(partner_exons)),
               rep(y_ctrl - exon_h / 2, nrow(partner_exons))),
      ymax = c(rep(y_sample + exon_h / 2, nrow(partner_exons)),
               rep(y_ctrl + exon_h / 2, nrow(partner_exons))),
      exon_type = rep(ifelse(partner_exons$is_focal, "focal", "other"), 2),
      stringsAsFactors = FALSE
    )
  )

  # sample plot data
  sample_plot <- merge(
    sample_df,
    partner_exons[, c("partner_pos", "xmid", "is_focal", "exon_start", "exon_end")],
    by.x = partner_col,
    by.y = "partner_pos",
    all.x = TRUE,
    sort = FALSE
  )

  sample_plot$x <- x_anchor + exon_w
  sample_plot$xend <- vapply(
    seq_len(nrow(sample_plot)),
    function(i) {
      pos_i <- sample_plot[[partner_col]][i]
      xm_i <- sample_plot$xmid[i]
      es_i <- sample_plot$exon_start[i]
      ee_i <- sample_plot$exon_end[i]

      map_pos_to_exon_x(
        pos = pos_i,
        xm = xm_i,
        exon_start = es_i,
        exon_end = ee_i,
        exon_w = exon_w
      )
    },
    numeric(1)
  )
  sample_plot$y <- y_sample + exon_h / 2
  sample_plot$yend <- y_sample + exon_h / 2
  sample_plot$label <- as.character(sample_plot$junction_count)

  # 去掉零长度/非法弧线
  sample_plot <- sample_plot[
    is.finite(sample_plot$x) &
      is.finite(sample_plot$xend) &
      abs(sample_plot$x - sample_plot$xend) > 1e-8,
    ,
    drop = FALSE
  ]

  # controls plot data
  others_plot <- merge(
    others_sum,
    partner_exons[, c("partner_pos", "xmid", "is_focal", "exon_start", "exon_end")],
    by.x = "partner",
    by.y = "partner_pos",
    all.x = TRUE,
    sort = FALSE
  )

  others_plot$x <- x_anchor + exon_w
  others_plot$xend <- vapply(
    seq_len(nrow(others_plot)),
    function(i) {
      map_pos_to_exon_x(
        pos = others_plot$partner[i],
        xm = others_plot$xmid[i],
        exon_start = others_plot$exon_start[i],
        exon_end = others_plot$exon_end[i],
        exon_w = exon_w
      )
    },
    numeric(1)
  )
  others_plot$y <- y_ctrl + exon_h / 2
  others_plot$yend <- y_ctrl + exon_h / 2
  others_plot$label <- paste0(round(others_plot$mean), " ± ", round(others_plot$sd))

  others_plot <- others_plot[
    is.finite(others_plot$x) &
      is.finite(others_plot$xend) &
      abs(others_plot$x - others_plot$xend) > 1e-8,
    ,
    drop = FALSE
  ]

  # labels
  if (nrow(sample_plot) > 0) {
    sample_plot <- sample_plot[order(sample_plot$xend), , drop = FALSE]
    sample_plot$label_x <- (sample_plot$x + sample_plot$xend) / 2
    sample_plot$label_y <- y_sample + seq(0.32, by = 0.22, length.out = nrow(sample_plot))
  }

  if (nrow(others_plot) > 0) {
    others_plot <- others_plot[order(others_plot$xend), , drop = FALSE]
    others_plot$label_x <- (others_plot$x + others_plot$xend) / 2
    others_plot$label_y <- y_ctrl + seq(0.32, by = 0.22, length.out = nrow(others_plot))
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rect_df,
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = exon_type
      ),
      color = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_curve(
      data = sample_plot,
      ggplot2::aes(
        x = x, y = y,
        xend = xend, yend = yend,
        color = is_focal
      ),
      curvature = 0.25,
      linewidth = 0.7
    ) +
    ggplot2::geom_curve(
      data = others_plot,
      ggplot2::aes(
        x = x, y = y,
        xend = xend, yend = yend,
        color = is_focal
      ),
      curvature = 0.25,
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = sample_plot,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      size = 4.8
    ) +
    ggplot2::geom_text(
      data = others_plot,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      size = 4.8
    ) +
    ggplot2::annotate(
      "text",
      x = x_anchor - 1.2,
      y = y_sample + 0.05,
      label = sampleID,
      hjust = 1,
      size = 5.5
    ) +
    ggplot2::annotate(
      "text",
      x = x_anchor - 1.2,
      y = y_ctrl + 0.05,
      label = "Controls (mean ± SD)",
      hjust = 1,
      size = 5.5
    )

  if (isTRUE(showReference)) {
    ref_rect <- data.frame(
      xmin = tx_ref$xmid - exon_w,
      xmax = tx_ref$xmid + exon_w,
      ymin = y_ref - exon_h / 2,
      ymax = y_ref + exon_h / 2,
      is_anchor = seq_len(nrow(tx_ref)) == anchor_local_idx,
      stringsAsFactors = FALSE
    )

    p <- p +
      ggplot2::geom_rect(
        data = ref_rect,
        ggplot2::aes(
          xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax,
          fill = is_anchor
        ),
        color = "black",
        linewidth = 0.4,
        show.legend = FALSE
      ) +
      ggplot2::annotate(
        "text",
        x = x_anchor - 1.2,
        y = y_ref + 0.05,
        label = "Reference transcript",
        hjust = 1,
        size = 5.2
      )
  }

  p +
    ggplot2::scale_fill_manual(
      values = c(
        anchor = "#2B1E1A",
        focal = "#FFFFFF",
        other = "#BDBDBD",
        `TRUE` = "#2B1E1A",
        `FALSE` = "#222222"
      ),
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c(
        `TRUE` = "#7A5230",
        `FALSE` = "#7A5230"
      ),
      guide = "none"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(x_anchor - 1.8, max(tx_ref$xmid) + 1.4),
      ylim = c(if (isTRUE(showReference)) -1.8 else -0.5, 3.0),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(20, 30, 20, 90)
    )
}
