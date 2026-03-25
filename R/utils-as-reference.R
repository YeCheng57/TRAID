buildReferenceJunctions <- function(
    transcript_exons,
    anchor_exon,
    partner_exons
) {
  # 排序 transcript exon
  tx <- transcript_exons[order(transcript_exons$start), , drop = FALSE]

  # 找 anchor exon 在 transcript 中的位置
  idx_anchor <- which(
    tx$start == anchor_exon$start &
      tx$end == anchor_exon$end
  )

  if (length(idx_anchor) == 0) {
    idx_anchor <- which.min(abs(tx$start - anchor_exon$start))
  }

  # 找 partner exon 在 transcript 中的位置
  partner_idx <- sapply(seq_len(nrow(partner_exons)), function(i) {
    hit <- which(
      tx$start == partner_exons$start[i] &
        tx$end == partner_exons$end[i]
    )
    if (length(hit) == 0) {
      which.min(abs(tx$start - partner_exons$start[i]))
    } else {
      hit[1]
    }
  })

  # 需要覆盖的 exon 区间
  idx_min <- min(c(idx_anchor, partner_idx))
  idx_max <- max(c(idx_anchor, partner_idx))

  # canonical junction：相邻 exon
  ref_pairs <- data.frame(
    idx_from = seq(idx_min, idx_max - 1),
    idx_to = seq(idx_min + 1, idx_max),
    stringsAsFactors = FALSE
  )

  ref_pairs
}
