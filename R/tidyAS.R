#' Return tidy AS result table
#'
#' @param x An `ASResult` object.
#' @param site One of "start", "end", "both".
#'
#' @return A data.frame if one site is requested, or a list if `site = "both"`.
#' @export
tidyAS <- function(x, site = c("start", "end", "both")) {
  if (!inherits(x, "ASResult")) {
    stop("`x` must be an ASResult object.")
  }

  site <- match.arg(site)

  if (site == "start") {
    return(x$result_start)
  }
  if (site == "end") {
    return(x$result_end)
  }

  list(
    start = x$result_start,
    end = x$result_end
  )
}
