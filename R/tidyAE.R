#' Return tidy AE result table
#'
#' @param x An `AEResult` object.
#'
#' @return A data.frame.
#' @export
tidyAE <- function(x) {
  if (!inherits(x, "AEResult")) {
    stop("`x` must be an AEResult object.")
  }
  x$result
}
