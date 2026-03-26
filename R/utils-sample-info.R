checkSampleInfo <- function(sample_info) {
  if (is.null(sample_info)) {
    return(invisible(TRUE))
  }

  if (!is.data.frame(sample_info)) {
    stop("`sample_info` must be a data.frame.")
  }

  req <- c("sample_id", "sex", "Gene_diagnostic")
  miss <- setdiff(req, colnames(sample_info))
  if (length(miss) > 0) {
    stop("`sample_info` is missing required columns: ",
         paste(miss, collapse = ", "))
  }

  if (anyDuplicated(sample_info$sample_id)) {
    stop("`sample_info$sample_id` must be unique.")
  }

  invisible(TRUE)
}

validateSamplesAgainstInfo <- function(samples, sample_info, input_name = "input") {
  if (is.null(sample_info)) {
    return(invisible(TRUE))
  }

  checkSampleInfo(sample_info)

  samples <- unique(as.character(samples))
  missing <- setdiff(samples, as.character(sample_info$sample_id))

  if (length(missing) > 0) {
    stop(
      "The following samples in ", input_name,
      " are not found in `sample_info$sample_id`: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}
