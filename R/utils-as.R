.compute_junction_psi <- function(df, site = c("start", "end")) {
  site <- match.arg(site)

  group_col <- if (site == "start") "start_id" else "end_id"

  total_df <- stats::aggregate(
    df$count,
    by = list(
      sample_id = df$sample_id,
      site_id = df[[group_col]]
    ),
    FUN = sum
  )
  colnames(total_df)[colnames(total_df) == "x"] <- "total_site_count"

  out <- merge(
    df,
    total_df,
    by.x = c("sample_id", group_col),
    by.y = c("sample_id", "site_id"),
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
    d$n_others_non_missing <- 0L
    d$reason_na <- NA_character_

    if (n < 2) {
      d$reason_na <- "insufficient_others"
      return(d)
    }

    for (i in seq_len(n)) {
      others <- d$psi[-i]
      others <- others[!is.na(others)]

      d$n_others_non_missing[i] <- length(others)

      if (length(others) < 2) {
        d$reason_na[i] <- "insufficient_others"
        next
      }

      mu <- mean(others)
      sigma <- stats::sd(others)

      d$psi_mean_others[i] <- mu
      d$psi_sd_others[i] <- sigma
      d$delta_psi[i] <- d$psi[i] - mu

      if (is.na(sigma) || sigma == 0) {
        d$reason_na[i] <- "zero_sd"
        next
      }

      d$zscore[i] <- d$delta_psi[i] / sigma
      d$pvalue[i] <- 2 * stats::pnorm(-abs(d$zscore[i]))
      d$reason_na[i] <- NA_character_
    }

    d
  })

  out <- do.call(rbind, res_list)
  out$padj <- stats::p.adjust(out$pvalue, method = pAdjustMethod)
  rownames(out) <- NULL
  out
}
