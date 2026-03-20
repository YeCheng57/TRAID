.compute_junction_psi <- function(df, site = c("donor", "acceptor")) {
  site <- match.arg(site)

  group_col <- if (site == "donor") "donor_id" else "acceptor_id"

  total_df <- stats::aggregate(
    count ~ sample_id + .data[[group_col]],
    data = df,
    FUN = sum
  )
  colnames(total_df)[colnames(total_df) == "count"] <- "total_site_count"

  out <- merge(
    df,
    total_df,
    by.x = c("sample_id", group_col),
    by.y = c("sample_id", group_col),
    all.x = TRUE,
    sort = FALSE
  )

  out$junction_count <- out$count
  out$psi <- out$junction_count / out$total_site_count
  out$site_type <- site
  out$site_id <- out[[group_col]]

  out
}


.compute_as_statistics_leave_one_out <- function(df, pAdjustMethod = "bonferroni") {
  split_idx <- split(seq_len(nrow(df)), df$junction_id)

  res_list <- lapply(split_idx, function(idx) {
    d <- df[idx, , drop = FALSE]

    n <- nrow(d)
    d$psi_mean_others <- NA_real_
    d$psi_sd_others <- NA_real_
    d$delta_psi <- NA_real_
    d$zscore <- NA_real_
    d$pvalue <- NA_real_

    if (n < 2) {
      return(d)
    }

    for (i in seq_len(n)) {
      others <- d$psi[-i]
      others <- others[!is.na(others)]

      if (length(others) < 2) {
        next
      }

      mu <- mean(others)
      sigma <- stats::sd(others)

      d$psi_mean_others[i] <- mu
      d$psi_sd_others[i] <- sigma
      d$delta_psi[i] <- d$psi[i] - mu

      if (!is.na(sigma) && sigma > 0) {
        d$zscore[i] <- d$delta_psi[i] / sigma
        d$pvalue[i] <- 2 * stats::pnorm(-abs(d$zscore[i]))
      }
    }

    d
  })

  out <- do.call(rbind, res_list)
  out$padj <- stats::p.adjust(out$pvalue, method = pAdjustMethod)
  rownames(out) <- NULL
  out
}
