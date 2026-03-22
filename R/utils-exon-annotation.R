.resolve_exon_annotation <- function(
    exon_annotation = NULL,
    gtf = NULL,
    cache = TRUE
) {
  if (!is.null(exon_annotation)) {
    .check_exon_annotation(exon_annotation)
    return(exon_annotation)
  }

  if (is.null(gtf)) {
    stop("Please provide either `exon_annotation` or `gtf`.")
  }

  if (!is.character(gtf) || length(gtf) != 1L) {
    stop("`gtf` must be a character scalar.")
  }

  is_url <- grepl("^https?://", gtf)

  cache_file <- NULL
  if (isTRUE(cache)) {
    key <- gsub("[^A-Za-z0-9]", "_", gtf)
    cache_file <- file.path(tempdir(), paste0("TRAID_exons_", key, ".rds"))

    if (file.exists(cache_file)) {
      exons <- readRDS(cache_file)
      .check_exon_annotation(exons)
      return(exons)
    }
  }

  gtf_file <- gtf

  if (is_url) {
    message("Downloading GTF...")
    ext <- if (grepl("\\.gz$", gtf)) ".gtf.gz" else ".gtf"
    gtf_file <- tempfile(fileext = ext)
    utils::download.file(gtf, destfile = gtf_file, mode = "wb", quiet = TRUE)
  }

  exons <- buildExonAnnotation(gtf_file)

  if (isTRUE(cache) && !is.null(cache_file)) {
    saveRDS(exons, cache_file)
  }

  exons
}
