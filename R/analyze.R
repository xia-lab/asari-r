# Corresponds to Python asari/analyze.py: single file mass track statistics and automatic peak height estimation.

# Call the database adapter method to avoid analyze.R's reverse dependence on the internal auxiliary function of annotate.R.
.analyze_call <- function(object, method, ...) {
  callback <- if (is.environment(object)) object[[method]] else object[[method]]
  if (!is.function(callback)) stop("Analysis adapter lacks method: ", method)
  callback(...)
}

# Read MS level and polarity from mzML or test injection table.
.analyze_scan_header <- function(infile, parameters) {
  if (!is.null(parameters$scan_header)) return(parameters$scan_header)
  if (!requireNamespace("mzR", quietly = TRUE)) stop("Analyzing mzML requires mzR.")
  handle <- mzR::openMSfile(infile)
  on.exit(mzR::close(handle), add = TRUE)
  mzR::header(handle)
}

# Corresponds to analyze_single_sample: prints statistics and estimates mass deviations using a known database.
analyze_single_sample <- function(infile, parameters = list()) {
  cat("Analysis of ", infile, "\n\n", sep = "")
  stats <- get_file_masstrack_stats(infile, parameters, return_sample = FALSE)
  factory <- getOption("asariR.ext_Experiment_factory")
  if (!is.function(factory)) factory <- get0("ext_Experiment", mode = "function", inherits = TRUE)
  if (!is.function(factory)) stop("ext_Experiment is unavailable.")
  experiment <- factory(list(), parameters)
  experiment$load_annotation_db()
  ratio <- .analyze_call(
    experiment$KCD, "evaluate_mass_accuracy_ratio",
    stats[[1L]], stats[[2L]], mz_tolerance_ppm = 10
  )
  cat(sprintf("  Mass accuracy is estimated as %2.1f ppm.\n\n", ratio * 1e6))
  invisible(NULL)
}

# Corresponds to get_file_masstrack_stats: extract the track, isotope anchor point and print the file summary.
get_file_masstrack_stats <- function(infile, parameters, return_sample = FALSE) {
  header <- .analyze_scan_header(infile, parameters)
  levels <- as.integer(header$msLevel)
  ms1 <- sum(levels == 1L)
  ms2 <- sum(levels == 2L)
  positive <- if ("polarity" %in% names(header)) {
    sum(levels == 1L & header$polarity == 1, na.rm = TRUE)
  } else 0L
  ms1_polarity <- if ("polarity" %in% names(header)) header$polarity[levels == 1L] else rep(0L, ms1)
  ionization_mode <- if (length(ms1_polarity) && tail(ms1_polarity, 1L) == 1) "pos" else "neg"

  extractor <- parameters$extract_mass_tracks
  if (!is.function(extractor)) extractor <- get0("extract_massTracks_", mode = "function", inherits = TRUE)
  if (!is.function(extractor)) extractor <- get0("extract_mass_tracks", mode = "function", inherits = TRUE)
  if (!is.function(extractor)) stop("Mass-track extractor is unavailable.")
  extracted <- extractor(
    infile, mz_tolerance_ppm = parameters$mz_tolerance_ppm,
    min_intensity = parameters$min_intensity_threshold,
    min_timepoints = parameters$min_timepoints,
    min_peak_height = parameters$min_peak_height
  )
  tracks <- lapply(seq_along(extracted$tracks), function(ii) {
    list(id_number = ii - 1L, mz = extracted$tracks[[ii]][[1L]], intensity = extracted$tracks[[ii]][[2L]])
  })
  sample <- list(
    sample_id = infile, input_file = infile, ion_mode = ionization_mode,
    list_scan_numbers = extracted$rt_numbers,
    list_retention_time = extracted$rt_times,
    list_mass_tracks = tracks
  )
  max_scan_number <- length(sample$list_scan_numbers)
  pair_finder <- parameters$find_mzdiff_pairs
  if (!is.function(pair_finder)) pair_finder <- get0(".constructors_find_mzdiff_pairs", mode = "function", inherits = TRUE)
  if (!is.function(pair_finder)) stop("m/z difference pair finder is unavailable.")
  pairs <- pair_finder(tracks, list_mz_diff = c(1.003355), mz_tolerance_ppm = parameters$mz_tolerance_ppm)
  pairs <- match_mzdiff_pairs_by_rt(
    pairs, tracks, max_scan_number, rt_window_perc = 0.05, min_scans_window = 10
  )
  landmarks <- flatten_tuplelist(pairs)
  all_mz <- vapply(tracks, `[[`, 0, "mz")
  landmark_tracks <- lapply(landmarks, function(index0) tracks[[as.integer(index0) + 1L]])
  heights <- vapply(landmark_tracks, function(track) max(track$intensity), 0)
  # Preserve NA summary without isotope anchor; automatic peak height path will give explicit error later.
  max_height <- if (length(heights)) as.integer(max(heights)) else NA_integer_
  min_height <- if (length(heights)) as.integer(min(heights)) else NA_integer_
  sample$anchor_mz_pairs <- pairs
  sample$number_anchor_mz_pairs <- length(pairs)
  sample$landmark_min_height <- min_height
  sample$landmark_max_height <- max_height

  cat("Total number of MS1 spectra: ", max_scan_number, "\n", sep = "")
  cat("of which ", positive, if (positive > 1L) " are" else " is", " positive ionization mode.\n", sep = "")
  cat("Assuming ionization mode is ", ionization_mode, ".\n\n", sep = "")
  cat("Total number of MS2 spectra: ", ms2, "\n", sep = "")
  max_rt <- if (length(sample$list_retention_time)) max(sample$list_retention_time) else NA_real_
  cat(sprintf("Maximum retention time (sec): %f\n", max_rt))
  if (length(all_mz)) {
    cat(sprintf(
      "m/z range: (min %f, median %f, max %f)\n\n",
      min(all_mz), stats::median(all_mz), max(all_mz)
    ))
  } else {
    cat("m/z range: no mass tracks were extracted.\n\n")
  }
  cat("Found ", length(tracks), " mass tracks.\n", sep = "")
  cat("Found ", length(pairs), " 12C/13C isotopic pairs as landmarks.\n", sep = "")
  cat("Max intensity in any landmark track:  ", format(max_height, big.mark = ",", scientific = FALSE), "\n", sep = "")
  cat("Minimal height of landmark tracks:  ", format(min_height, big.mark = ",", scientific = FALSE), " \n\n", sep = "")
  if (isTRUE(return_sample)) return(sample)
  if (!length(heights)) {
    stop("No isotopic landmark pair was found for peak-height estimation.", call. = FALSE)
  }
  list(landmarks, ionization_mode, min_height)
}

# Corresponds to match_mzdiff_pairs_by_rt: The difference between the two apex scans is strictly smaller than the dynamic window and is retained.
match_mzdiff_pairs_by_rt <- function(
    matched_mz_pairs, list_mass_tracks, max_scan_number,
    rt_window_perc = 0.05, min_scans_window = 10) {
  window <- max(min_scans_window, max_scan_number * rt_window_perc)
  tracks <- stats::setNames(
    lapply(list_mass_tracks, `[[`, "intensity"),
    vapply(list_mass_tracks, function(track) as.character(track$id_number), "")
  )
  Filter(function(pair) {
    apex1 <- which.max(tracks[[as.character(pair[[1L]])]]) - 1L
    apex2 <- which.max(tracks[[as.character(pair[[2L]])]]) - 1L
    abs(apex1 - apex2) < window
  }, matched_mz_pairs)
}

# Corresponds to __wrapped_get_file_masstrack_stats: captures all single file exceptions.
`__wrapped_get_file_masstrack_stats` <- function(job) {
  infile <- job[[1L]]
  tryCatch(
    get_file_masstrack_stats(infile, job[[2L]], return_sample = FALSE),
    error = function(error) {
      cat("Error in analyzing ", infile, "\n", sep = "")
      NULL
    }
  )
}

# Corresponds to estimate_min_peak_height: sample the file and take the median half of the lowest peak height of the landmark.
estimate_min_peak_height <- function(list_input_files, parameters) {
  lower_bound <- max(
    parameters$min_min_peak_height,
    parameters$min_peak_height / parameters$dynamic_range
  )
  number <- parameters$num_files_to_check
  selected <- if (is.null(number)) {
    list_input_files
  } else if (length(list_input_files) <= number) {
    list_input_files
  } else if (is.integer(number) || (is.numeric(number) && number == as.integer(number))) {
    sample(list_input_files, as.integer(number))
  } else if (is.numeric(number) && number >= 0 && number <= 1) {
    to_select <- max(as.integer(number * length(list_input_files)), 10L)
    sample(list_input_files, min(to_select, length(list_input_files)))
  } else {
    stop("num_files_to_use should be int or float between 0 and 1 or None")
  }
  cat("Estimating parameter for min peak_height based on ", paste(selected, collapse = ", "), "\n")
  jobs <- lapply(selected, function(infile) list(infile, parameters))
  results <- bulk_process(
    `__wrapped_get_file_masstrack_stats`, jobs,
    jobs_per_worker = parameters$multicores
  )
  estimated <- numeric()
  for (ii in seq_along(results)) {
    if (!is.null(results[[ii]])) estimated <- c(estimated, results[[ii]][[3L]])
    else cat("Error in analyzing ", selected[[ii]], "\n", sep = "")
  }
  if (length(estimated) == 0L) stop("No input file produced a peak-height estimate.")
  recommended <- as.integer(0.5 * stats::median(estimated))
  cat("Estimated parameter for min peak_height is ", recommended, " \n\n", sep = "")
  max(recommended, lower_bound)
}
