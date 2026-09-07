# Corresponds to Python asari/tools/extract_masstracks.py: the old version of the independent mass track export tool.

# Corresponds to read_project_dir: no sorting, retaining directory enumeration order and filtering by substring.
read_project_dir_masstracks <- function(directory, file_pattern = ".mzML") {
  entries <- list.files(directory, full.names = FALSE, all.files = FALSE, no.. = TRUE)
  file.path(directory, entries[grepl(file_pattern, entries, fixed = TRUE)])
}

# Provide an explicit alias of the original function name; after workflow.R is loaded, it may overwrite the main process function with the same name.
read_project_dir_extract_masstracks <- read_project_dir_masstracks

# Read fields from the environment or list.
.extract_masstracks_get <- function(object, name) {
  if (is.environment(object)) return(object[[name]])
  object[[name]]
}

# Corresponds to process_single_file: extract mass tracks and export summary statistics for each track.
process_single_file <- function(infile, outdir = "") {
  # The original Python uses the outdated SimpleSample(input_file=...) interface; the R version allows the injection of factories with the same function.
  factory <- getOption("asariR.simple_sample_from_file")
  if (!is.function(factory)) {
    factory <- get0("SimpleSample_from_file", mode = "function", inherits = TRUE)
  }
  if (!is.function(factory)) {
    stop(
      "process_single_file requires options(asariR.simple_sample_from_file = function(path) ...).",
      call. = FALSE
    )
  }
  sample <- factory(infile)
  extractor <- .extract_masstracks_get(sample, "get_mass_tracks_")
  if (!is.function(extractor)) stop("Sample object lacks get_mass_tracks_.")
  extractor(mz_tolerance_ppm = 5, min_intensity = 100, min_timepoints = 5)
  tracks <- .extract_masstracks_get(sample, "list_mass_tracks")
  rows <- lapply(tracks, function(track) {
    intensity <- as.numeric(track$intensity)
    positive <- intensity[intensity > 0]
    median_positive <- if (length(positive)) as.integer(stats::median(positive)) else NA_integer_
    c(
      as.character(track$id_number), as.character(track$mz),
      as.character(length(intensity)), as.character(max(intensity)),
      as.character(median_positive)
    )
  })
  lines <- c(
    "id_number\tmz\tnumber_scans\tmax_intensity\tmedian_intensity",
    vapply(rows, paste, collapse = "\t", FUN.VALUE = "")
  )
  outfile <- file.path(outdir, sub(".mzML", "_massTrack.tsv", basename(infile), fixed = TRUE))
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}
