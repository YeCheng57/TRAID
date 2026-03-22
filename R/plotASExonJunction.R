#' Plot exon-aware junction cartoon for one AS event
#'
#' @param x An `ASResult` object, preferably generated with `returnAll = TRUE`.
#' @param sampleID Character scalar specifying the target sample.
#' @param chrom Chromosome of the focal junction.
#' @param start Junction start coordinate.
#' @param end Junction end coordinate.
#' @param site One of "start" or "end".
#' @param exon_annotation A data.frame containing exon-level annotation.
#' It must include the following columns:
#' \describe{
#'   \item{chrom}{Character. Chromosome name (e.g. `"1"` or `"chr1"`).}
#'   \item{start}{Integer. Exon start coordinate (1-based, inclusive).}
#'   \item{end}{Integer. Exon end coordinate (1-based, inclusive).}
#'   \item{strand}{Character. Strand information (`"+"`, `"-"`, or `"*"`).}
#'   \item{gene_id}{Character. Gene identifier.}
#'   \item{gene_name}{Character. Gene symbol or gene name.}
#'   \item{transcript_id}{Character. Transcript identifier.}
#'   \item{exon_id}{Character. Unique exon identifier.}
#'   \item{exon_number}{Integer. Exon number within transcript (optional but recommended).}
#' }
#' @param gtf Optional local path or URL to a GTF/GTF.GZ file. Used only if
#' `exon_annotation` is not provided.
#' @param cache Logical; whether to cache parsed exon annotation when `gtf` is used.
#' @param topN Integer; number of competing junctions to display.
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
    topN = 3
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

  gene_name <- NA_character_
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
    gene_name = if (is.na(gene_name)) NULL else gene_name,
    exon_annotation = exons
  )

  tx_exons <- exon_ctx$transcript_exons
  if (nrow(tx_exons) == 0) {
    stop("No transcript exons found for plotting.")
  }

  sample_df <- event_df[event_df$sample_id == sampleID, , drop = FALSE]
  others_df <- event_df[event_df$sample_id != sampleID, , drop = FALSE]

  if (nrow(sample_df) == 0) {
    stop("Target sample has no junctions at this focal site.")
  }

  if (nrow(others_df) > 0) {
    mean_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = others_df[[partner_col]]),
      FUN = mean,
      na.rm = TRUE
    )
    sd_df <- stats::aggregate(
      others_df$junction_count,
      by = list(partner = others_df[[partner_col]]),
      FUN = stats::sd,
      na.rm = TRUE
    )
    colnames(mean_df)[2] <- "mean_count"
    colnames(sd_df)[2] <- "sd_count"
    others_sum <- merge(mean_df, sd_df, by = "partner", all = TRUE)
  } else {
    others_sum <- data.frame(
      partner = numeric(0),
      mean_count = numeric(0),
      sd_count = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  sample_rank <- sample_df[order(sample_df$junction_count, decreasing = TRUE), , drop = FALSE]
  keep_partners <- unique(c(focal_partner, head(sample_rank[[partner_col]], topN)))

  sample_df <- sample_df[sample_df[[partner_col]] %in% keep_partners, , drop = FALSE]
  others_sum <- others_sum[others_sum$partner %in% keep_partners, , drop = FALSE]

  tx_exons <- tx_exons[order(tx_exons$start, tx_exons$end), , drop = FALSE]
  tx_exons$xmid <- seq_len(nrow(tx_exons)) * 2.4

  locate_exon_for_pos <- function(pos, ex_df) {
    hit <- ex_df[ex_df$start <= pos & ex_df$end >= pos, , drop = FALSE]
    if (nrow(hit) > 0) {
      return(hit[1, , drop = FALSE])
    }
    mid <- (ex_df$start + ex_df$end) / 2
    ex_df[which.min(abs(mid - pos)), , drop = FALSE]
  }

  anchor_exon <- locate_exon_for_pos(anchor_pos, tx_exons)

  partner_exons <- lapply(keep_partners, function(pos) {
    locate_exon_for_pos(pos, tx_exons)
  })
  partner_exons <- do.call(rbind, partner_exons)
  rownames(partner_exons) <- NULL
  partner_exons$partner_pos <- keep_partners
  partner_exons$is_focal <- partner_exons$partner_pos == focal_partner

  exon_w <- 0.55
  exon_h <- 0.20
  y_sample <- 2
  y_ctrl <- 0

  anchor_rect <- data.frame(
    panel = c("Sample", "Controls"),
    xmin = anchor_exon$xmid - exon_w,
    xmax = anchor_exon$xmid + exon_w,
    ymin = c(y_sample - exon_h / 2, y_ctrl - exon_h / 2),
    ymax = c(y_sample + exon_h / 2, y_ctrl + exon_h / 2),
    exon_type = "anchor",
    partner_pos = NA_real_,
    xmid = anchor_exon$xmid,
    stringsAsFactors = FALSE
  )

  partner_rect <- data.frame(
    panel = rep(c("Sample", "Controls"), each = nrow(partner_exons)),
    xmin = rep(partner_exons$xmid - exon_w, 2),
    xmax = rep(partner_exons$xmid + exon_w, 2),
    ymin = c(
      rep(y_sample - exon_h / 2, nrow(partner_exons)),
      rep(y_ctrl - exon_h / 2, nrow(partner_exons))
    ),
    ymax = c(
      rep(y_sample + exon_h / 2, nrow(partner_exons)),
      rep(y_ctrl + exon_h / 2, nrow(partner_exons))
    ),
    exon_type = rep(ifelse(partner_exons$is_focal, "focal", "other"), 2),
    partner_pos = rep(partner_exons$partner_pos, 2),
    xmid = rep(partner_exons$xmid, 2),
    stringsAsFactors = FALSE
  )

  rect_df <- rbind(anchor_rect, partner_rect)

  sample_plot <- merge(
    sample_df,
    partner_exons[, c("partner_pos", "xmid", "is_focal")],
    by.x = partner_col,
    by.y = "partner_pos",
    all.x = TRUE,
    sort = FALSE
  )

  sample_plot$panel <- "Sample"
  sample_plot$x <- anchor_exon$xmid
  sample_plot$xend <- sample_plot$xmid
  sample_plot$y <- y_sample + exon_h / 2
  sample_plot$yend <- y_sample + exon_h / 2
  sample_plot$label <- as.character(sample_plot$junction_count)
  sample_plot$label_x <- (sample_plot$x + sample_plot$xend) / 2
  sample_plot$label_y <- y_sample + 0.45

  others_plot <- merge(
    others_sum,
    partner_exons[, c("partner_pos", "xmid", "is_focal")],
    by.x = "partner",
    by.y = "partner_pos",
    all.x = TRUE,
    sort = FALSE
  )

  others_plot$panel <- "Controls"
  others_plot$x <- anchor_exon$xmid
  others_plot$xend <- others_plot$xmid
  others_plot$y <- y_ctrl + exon_h / 2
  others_plot$yend <- y_ctrl + exon_h / 2
  others_plot$label <- paste0(
    round(others_plot$mean_count, 0),
    " ± ",
    round(others_plot$sd_count, 0)
  )
  others_plot$label_x <- (others_plot$x + others_plot$xend) / 2
  others_plot$label_y <- y_ctrl + 0.45

  rank_heights <- function(n, base_y) {
    base_y + seq(0.12, by = 0.12, length.out = n)
  }

  if (nrow(sample_plot) > 0) {
    sample_plot <- sample_plot[order(sample_plot$xend), , drop = FALSE]
    sample_plot$label_y <- rank_heights(nrow(sample_plot), y_sample + 0.30)
  }

  if (nrow(others_plot) > 0) {
    others_plot <- others_plot[order(others_plot$xend), , drop = FALSE]
    others_plot$label_y <- rank_heights(nrow(others_plot), y_ctrl + 0.30)
  }

  ggplot2::ggplot() +
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
        x = x + exon_w * 0.9,
        y = y,
        xend = xend - exon_w * 0.9,
        yend = yend,
        color = is_focal
      ),
      curvature = -0.25,
      linewidth = 0.7
    ) +
    ggplot2::geom_curve(
      data = others_plot,
      ggplot2::aes(
        x = x + exon_w * 0.9,
        y = y,
        xend = xend - exon_w * 0.9,
        yend = yend,
        color = is_focal
      ),
      curvature = -0.25,
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
      x = min(c(anchor_exon$xmid, partner_exons$xmid)) - 1.8,
      y = y_sample + 0.10,
      label = sampleID,
      hjust = 1,
      size = 5.5
    ) +
    ggplot2::annotate(
      "text",
      x = min(c(anchor_exon$xmid, partner_exons$xmid)) - 1.8,
      y = y_ctrl + 0.10,
      label = "Controls (mean ± SD)",
      hjust = 1,
      size = 5.5
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        anchor = "#222222",
        focal = "#FFFFFF",
        other = "#BDBDBD"
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
      xlim = c(
        min(c(anchor_exon$xmid, partner_exons$xmid)) - 2.5,
        max(c(anchor_exon$xmid, partner_exons$xmid)) + 1.2
      ),
      ylim = c(-0.5, 3.0),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(20, 30, 20, 90)
    )
}
