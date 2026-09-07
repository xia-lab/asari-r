# Additional public entry points expose the completed pure-R modules through stable, verifiable functions.

# Checks whether an argument is a single positive number and returns the numeric form.
.asari_public_positive_number <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
    stop(name, " must be one positive number.", call. = FALSE)
  }
  as.numeric(value)
}

# Read feature table from an explicit path or data.frame, without looking for any default file.
.asari_public_feature_table <- function(value, name = "feature_table") {
  if (is.data.frame(value)) return(value)
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop(name, " must be a data.frame or one explicit TSV file path.", call. = FALSE)
  }
  path <- path.expand(value)
  if (!file.exists(path) || dir.exists(path)) {
    stop(name, " does not exist or is not a file: ", path, call. = FALSE)
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Check the columns that feature table must have, and the missing names are listed directly in the error.
.asari_public_require_columns <- function(table, columns, name = "feature_table") {
  missing <- setdiff(columns, names(table))
  if (length(missing) > 0L) {
    stop(
      name, " is missing required column(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  invisible(TRUE)
}

# Prepare an explicit output file and create its parent directory.
.asari_public_output_file <- function(output, name = "output") {
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(trimws(output))) {
    stop(name, " must be one explicit file path.", call. = FALSE)
  }
  path <- path.expand(output)
  parent <- dirname(path)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", parent, call. = FALSE)
  }
  normalizePath(path, mustWork = FALSE)
}

# Check and return an explicit mzML file.
.asari_public_one_mzml <- function(input) {
  files <- .asari_collect_input_files(input, recursive = FALSE)
  if (length(files) != 1L) {
    stop("This function requires exactly one explicit .mzML file.", call. = FALSE)
  }
  files[[1L]]
}

#' Use the DIMS process to process direct injection mass spectrum data
#'
#' This function reuses the input check, mass track extraction and result export of [asari_process()],
#' But fixed use `DIMS` process and turn off retention time alignment.
#'
#' @inheritParams asari_process
#' @return An `asari_result` object. The main DIMS result is stored in `full_feature_table`.
#' @export
asari_process_dims <- function(
    input,
    output = "asari-output",
    project_name = "asari_dims",
    mode = NULL,
    ppm = NULL,
    multicores = NULL,
    database_mode = NULL,
    recursive = FALSE,
    parameters = list()) {
  # DIMS does not use chromatogram peak boundaries or cross-sample retention time alignment.
  parameters$workflow <- "DIMS"
  parameters$rt_align_on <- FALSE
  asari_process(
    input = input,
    output = output,
    project_name = project_name,
    mode = mode,
    ppm = ppm,
    multicores = multicores,
    rt_align = FALSE,
    database_mode = database_mode,
    recursive = recursive,
    parameters = parameters
  )
}

#' Analyze the scans and mass tracks in an mzML file
#'
#' This function only performs local mzML reading, mass track extraction and isotope anchor point statistics.
#' Known compound databases are not connected and external annotation is not performed.
#'
#' @param input An explicit `.mzML` file.
#' @param ppm m/z tolerance in ppm.
#' @param min_intensity Minimum data point intensity.
#' @param min_timepoints The minimum number of data points contained in a mass track.
#' @param min_peak_height Minimum peak height.
#' @param parameters A named override list of other advanced parameters.
#' @return Contains scan summary, anchor statistics, mass track table, and named list of raw trajectories.
#' @export
asari_analyze <- function(
    input,
    ppm = 5,
    min_intensity = 1000,
    min_timepoints = 6L,
    min_peak_height = 100000,
    parameters = list()) {
  # All analyzes start from the unique file given by the user and do not read other mzML in the working directory.
  infile <- .asari_public_one_mzml(input)
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Analyzing mzML requires the Bioconductor package mzR.", call. = FALSE)
  }
  effective <- .asari_merge_parameters(parameters)
  effective$mz_tolerance_ppm <- .asari_public_positive_number(ppm, "ppm")
  effective$min_intensity_threshold <- .asari_public_positive_number(
    min_intensity, "min_intensity"
  )
  effective$min_timepoints <- as.integer(
    .asari_public_positive_number(min_timepoints, "min_timepoints")
  )
  effective$min_peak_height <- .asari_public_positive_number(
    min_peak_height, "min_peak_height"
  )

  # Read the scan header only once and inject it into the existing analysis module to avoid opening large files repeatedly.
  header <- .analyze_scan_header(infile, effective)
  effective$scan_header <- header
  sample <- get_file_masstrack_stats(infile, effective, return_sample = TRUE)
  tracks <- sample$list_mass_tracks
  track_table <- if (length(tracks) == 0L) {
    data.frame(
      id_number = integer(), mz = numeric(), number_scans = integer(),
      max_intensity = numeric(), median_positive_intensity = numeric()
    )
  } else {
    do.call(rbind, lapply(tracks, function(track) {
      intensity <- as.numeric(track$intensity)
      positive <- intensity[intensity > 0]
      data.frame(
        id_number = as.integer(track$id_number),
        mz = as.numeric(track$mz),
        number_scans = length(intensity),
        max_intensity = max(intensity),
        median_positive_intensity = if (length(positive)) {
          stats::median(positive)
        } else {
          NA_real_
        }
      )
    }))
  }

  # Returns structured results that are both suitable for quick viewing and sufficient to continue downstream calculations.
  levels <- as.integer(header$msLevel)
  result <- list(
    input_file = infile,
    ion_mode = sample$ion_mode,
    ms1_scans = sum(levels == 1L, na.rm = TRUE),
    ms2_scans = sum(levels == 2L, na.rm = TRUE),
    max_retention_time = if (length(sample$list_retention_time)) {
      max(sample$list_retention_time)
    } else {
      NA_real_
    },
    mass_track_count = length(tracks),
    isotope_pair_count = length(sample$anchor_mz_pairs),
    landmark_min_height = sample$landmark_min_height,
    landmark_max_height = sample$landmark_max_height,
    mass_track_table = track_table,
    mass_tracks = tracks,
    parameters = effective[names(effective) != "scan_header"]
  )
  class(result) <- "asari_analysis"
  result
}

#' Export mass track summary of an mzML file
#'
#' @inheritParams asari_analyze
#' @param output Optional explicit TSV output path; only data.frame is returned when it is `NULL`.
#' @return mass track summary data.frame.
#' @export
asari_extract_mass_tracks <- function(
    input,
    output = NULL,
    ppm = 5,
    min_intensity = 1000,
    min_timepoints = 6L,
    min_peak_height = 100000,
    parameters = list()) {
  # Reuse single-file analysis to ensure that only one set of trajectory extraction parameters and statistical criteria is implemented.
  analysis <- asari_analyze(
    input = input,
    ppm = ppm,
    min_intensity = min_intensity,
    min_timepoints = min_timepoints,
    min_peak_height = min_peak_height,
    parameters = parameters
  )
  table <- analysis$mass_track_table
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(
      table, outfile, sep = "\t", row.names = FALSE,
      col.names = TRUE, quote = FALSE, na = ""
    )
    attr(table, "output_file") <- outfile
  }
  attr(table, "analysis") <- analysis
  table
}

#' Read the project results that asariR has generated
#'
#' @param project A project directory, or `asari_result` returned by [asari_process()].
#' @param table To read the `"preferred"`, `"full"` or `"both"` table.
#' @param max_samples Optional maximum number of sample columns; `NULL` means reading all sample columns.
#' @return An `asari_project_results` object containing the project description and selected feature table.
#' @export
asari_read_results <- function(
    project,
    table = c("both", "preferred", "full"),
    max_samples = NULL) {
  table <- match.arg(table)
  project_dir <- if (inherits(project, "asari_result")) {
    project$project_dir
  } else {
    project
  }
  if (!is.character(project_dir) || length(project_dir) != 1L ||
      is.na(project_dir) || !dir.exists(path.expand(project_dir))) {
    stop("project must be an existing project directory or asari_result.", call. = FALSE)
  }
  project_dir <- normalizePath(path.expand(project_dir), mustWork = TRUE)
  description_file <- file.path(project_dir, "project.json")
  if (!file.exists(description_file)) {
    stop("Project is missing project.json: ", project_dir, call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Reading project.json requires jsonlite.", call. = FALSE)
  }
  if (!is.null(max_samples)) {
    if (!is.numeric(max_samples) || length(max_samples) != 1L ||
        is.na(max_samples) || max_samples < 0 || max_samples != as.integer(max_samples)) {
      stop("max_samples must be NULL or one non-negative integer.", call. = FALSE)
    }
    max_samples <- as.integer(max_samples)
  }

  # The first 11 columns of the asari table are fixed feature metadata, and the following columns are sample abundances.
  read_one <- function(path) {
    if (!file.exists(path)) stop("Project is missing feature table: ", path, call. = FALSE)
    value <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    if (!is.null(max_samples) && ncol(value) > 11L + max_samples) {
      value <- value[, seq_len(11L + max_samples), drop = FALSE]
    }
    value
  }
  project_description <- jsonlite::fromJSON(
    description_file, simplifyVector = FALSE
  )
  # Use the table name actually recorded in the project, compatible with the output_feature_table advanced override value.
  feature_table_name <- if (!is.null(project_description$output_feature_table)) {
    as.character(project_description$output_feature_table)
  } else {
    "Feature_table.tsv"
  }
  result <- list(
    project_dir = project_dir,
    project = project_description,
    preferred = NULL,
    full = NULL
  )
  if (table %in% c("both", "preferred")) {
    result$preferred <- read_one(file.path(
      project_dir, paste0("preferred_", feature_table_name)
    ))
  }
  if (table %in% c("both", "full")) {
    result$full <- read_one(file.path(
      project_dir, "export", paste0("full_", feature_table_name)
    ))
  }
  class(result) <- "asari_project_results"
  result
}

#' Compare two LC-MS feature tables
#'
#' The candidate with the smallest m/z deviation within the ppm and RT windows is selected in both directions, and only the best matches to each other are retained.
#'
#' @param left,right TSV path or data.frame.
#' @param ppm Positive m/z tolerance.
#' @param rt_tolerance Positive retention time tolerance, unit is seconds, using strict less than relationship.
#' @param left_rt_unit,right_rt_unit Each table has retention time unit.
#' @param verbose Whether to print the underlying two-way matching statistics.
#' @return One data.frame per row that best matches each other.
#' @export
asari_compare_features <- function(
    left,
    right,
    ppm = 5,
    rt_tolerance = 5,
    left_rt_unit = c("seconds", "minutes"),
    right_rt_unit = c("seconds", "minutes"),
    verbose = FALSE) {
  ppm <- .asari_public_positive_number(ppm, "ppm")
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  left_rt_unit <- match.arg(left_rt_unit)
  right_rt_unit <- match.arg(right_rt_unit)
  left_table <- .asari_public_feature_table(left, "left")
  right_table <- .asari_public_feature_table(right, "right")
  .asari_public_require_columns(left_table, c("mz", "rtime"), "left")
  .asari_public_require_columns(right_table, c("mz", "rtime"), "right")

  # Use a stable row number ID when there is no id_number to avoid requiring users to modify their own tables first.
  if (!("id_number" %in% names(left_table))) {
    left_table$id_number <- paste0("row", seq_len(nrow(left_table)))
  }
  if (!("id_number" %in% names(right_table))) {
    right_table$id_number <- paste0("row", seq_len(nrow(right_table)))
  }
  left_table$id_number <- as.character(left_table$id_number)
  right_table$id_number <- as.character(right_table$id_number)
  if (anyDuplicated(left_table$id_number) || anyDuplicated(right_table$id_number)) {
    stop("id_number values must be unique within each feature table.", call. = FALSE)
  }
  left_seconds <- as.numeric(left_table$rtime) * if (left_rt_unit == "minutes") 60 else 1
  right_seconds <- as.numeric(right_table$rtime) * if (right_rt_unit == "minutes") 60 else 1
  if (anyNA(c(as.numeric(left_table$mz), left_seconds)) ||
      anyNA(c(as.numeric(right_table$mz), right_seconds))) {
    stop("mz and rtime columns must contain numeric values without NA.", call. = FALSE)
  }
  left_features <- lapply(seq_len(nrow(left_table)), function(ii) list(
    id = left_table$id_number[[ii]], mz = as.numeric(left_table$mz[[ii]]),
    rtime = left_seconds[[ii]]
  ))
  right_features <- lapply(seq_len(nrow(right_table)), function(ii) list(
    id = right_table$id_number[[ii]], mz = as.numeric(right_table$mz[[ii]]),
    rtime = right_seconds[[ii]]
  ))

  # Directly call the two-way best matching algorithm that has been translated and tested def by def.
  matched <- if (isTRUE(verbose)) {
    bidirectional_best_match(left_features, right_features, ppm, rt_tolerance)
  } else {
    invisible(capture.output(
      value <- bidirectional_best_match(
        left_features, right_features, ppm, rt_tolerance
      )
    ))
    value
  }
  pairs <- matched[[1L]]
  if (length(pairs) == 0L) {
    result <- data.frame(
      left_id = character(), right_id = character(),
      left_mz = numeric(), right_mz = numeric(), ppm_error = numeric(),
      left_rtime = numeric(), right_rtime = numeric(), rt_delta = numeric()
    )
  } else {
    result <- do.call(rbind, lapply(pairs, function(pair) {
      li <- match(pair[[1L]], left_table$id_number)
      ri <- match(pair[[2L]], right_table$id_number)
      data.frame(
        left_id = pair[[1L]], right_id = pair[[2L]],
        left_mz = as.numeric(left_table$mz[[li]]),
        right_mz = as.numeric(right_table$mz[[ri]]),
        ppm_error = (as.numeric(right_table$mz[[ri]]) -
          as.numeric(left_table$mz[[li]])) /
          as.numeric(left_table$mz[[li]]) * 1e6,
        left_rtime = left_seconds[[li]], right_rtime = right_seconds[[ri]],
        rt_delta = right_seconds[[ri]] - left_seconds[[li]],
        stringsAsFactors = FALSE
      )
    }))
  }
  attr(result, "summary") <- c(
    left_features = nrow(left_table), right_features = nrow(right_table),
    mutual_matches = nrow(result), left_unmatched = nrow(left_table) - nrow(result),
    right_unmatched = nrow(right_table) - nrow(result)
  )
  result
}

#' Extract candidate peaks from feature table by target m/z
#'
#' @param feature_table asari feature table path or data.frame.
#' @param targets One or more target m/z.
#' @param ppm Positive m/z tolerance.
#' @param output Optional explicit TSV output path.
#' @return Hit feature rows; first column is query mass, with ppm error.
#' @export
asari_extract_targets <- function(feature_table, targets, ppm = 5, output = NULL) {
  table <- .asari_public_feature_table(feature_table)
  .asari_public_require_columns(table, "mz")
  ppm <- .asari_public_positive_number(ppm, "ppm")
  if (!is.numeric(targets) || length(targets) == 0L || anyNA(targets) ||
      any(targets <= 0)) {
    stop("targets must contain one or more positive numeric m/z values.", call. = FALSE)
  }
  feature_mz <- as.numeric(table$mz)
  if (anyNA(feature_mz)) stop("feature_table$mz must be numeric without NA.", call. = FALSE)
  mapping <- all_mass_paired_mapping(feature_mz, as.numeric(targets), ppm)
  pairs <- mapping$mapped
  if (length(pairs) == 0L) {
    result <- table[FALSE, , drop = FALSE]
    result <- cbind(
      query_target = numeric(), ppm_error = numeric(), result,
      stringsAsFactors = FALSE
    )
  } else {
    feature_index <- vapply(pairs, `[[`, integer(1), 1L)
    target_index <- vapply(pairs, `[[`, integer(1), 2L)
    query <- as.numeric(targets)[target_index]
    result <- cbind(
      query_target = query,
      ppm_error = (feature_mz[feature_index] - query) / query * 1e6,
      table[feature_index, , drop = FALSE],
      stringsAsFactors = FALSE
    )
  }
  attr(result, "unmatched_targets") <- as.numeric(targets)[mapping$list2_unmapped]
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(
      result, outfile, sep = "\t", row.names = FALSE,
      col.names = TRUE, quote = FALSE, na = ""
    )
    attr(result, "output_file") <- outfile
  }
  result
}

#' Generate a QC plot for a feature table
#'
#' @param feature_table asari full feature table path or data.frame.
#' @param output Explicit PDF output path.
#' @param height PDF height.
#' @param aspect PDF aspect ratio.
#' @return Generated PDF absolute path.
#' @export
asari_feature_qc <- function(feature_table, output, height = 12, aspect = 0.7) {
  table <- .asari_public_feature_table(feature_table)
  .asari_public_require_columns(
    table,
    c("snr", "peak_area", "cSelectivity", "goodness_fitting")
  )
  outfile <- .asari_public_output_file(output)
  asari_qc_plot(table, outfile = outfile, height = height, aspect = aspect)
  normalizePath(outfile, mustWork = TRUE)
}

#' Generate HTML QA report for an mzML file
#'
#' @param input An explicit `.mzML` file.
#' @param output Explicit HTML output path.
#' @param spikeins `NULL`, target list or target JSON file.
#' @return Generated HTML absolute path.
#' @export
asari_qc_report <- function(input, output, spikeins = NULL) {
  infile <- .asari_public_one_mzml(input)
  outfile <- .asari_public_output_file(output)
  generate_qc_report(list(infile, outfile, spikeins))
  normalizePath(outfile, mustWork = TRUE)
}

#' Generate MS1/MS2 scan summary PDF for an mzML file
#'
#' @param input An explicit `.mzML` file.
#' @param output Explicit PDF output path.
#' @param nspec_plot Number of sample scans to display.
#' @param offset_n Display starts from the vicinity of the scan number.
#' @param width,height PDF size.
#' @param title Figure title.
#' @return Generated PDF absolute path.
#' @export
asari_scan_summary <- function(
    input,
    output,
    nspec_plot = 50L,
    offset_n = 100L,
    width = 10,
    height = 6,
    title = "") {
  infile <- .asari_public_one_mzml(input)
  outfile <- .asari_public_output_file(output)
  summarize_ms_file(
    infile = infile, nspec_plot = nspec_plot, offset_n = offset_n,
    width = width, height = height, title = title, outfile = outfile
  )
  normalizePath(outfile, mustWork = TRUE)
}

#' View the available status of the current workflow
#'
#' @return A data.frame describing local direct use and external dependency status.
#' @export
asari_available_workflows <- function() {
  # These status descriptions document the public entry points available after installation.
  data.frame(
    workflow = c("LC", "DIMS", "GC", "LCMSMS", "LC_annotation", "dashboard", "RAW_conversion"),
    direct_use = rep(TRUE, 7L),
    public_function = c(
      "asari_process", "asari_process_dims", "asari_annotate_gc",
      "asari_extract_ms2 / asari_match_ms2 / asari_search_ms2",
      "asari_annotate_lc", "asari_dashboard", "asari_convert_raw"
    ),
    note = c(
      "Local core workflow",
      "Local direct-infusion workflow; use full_feature_table",
      "Requires explicit Kovats standards and an MSP/JSON EI library",
      "Local extraction and search; entropy uses a documented local fallback",
      "Mass/RT candidate annotation; hmdb4 first load needs Python pickle decoding",
      "Self-contained offline HTML; no Shiny or dashboard pickle required",
      "Requires ThermoRawFileParser and, on macOS/Linux, a compatible runtime"
    ),
    stringsAsFactors = FALSE
  )
}

# Prints a single-file analysis summary and does not automatically expand a potentially large mass-track list.
#' @export
print.asari_analysis <- function(x, ...) {
  invisible(list(...))
  cat("asariR single-file analysis.\n")
  cat("Input: ", x$input_file, "\n", sep = "")
  cat("Ion mode: ", x$ion_mode, "\n", sep = "")
  cat("MS1 scans: ", x$ms1_scans, "; MS2 scans: ", x$ms2_scans, "\n", sep = "")
  cat("Mass tracks: ", x$mass_track_count,
      "; isotope pairs: ", x$isotope_pair_count, "\n", sep = "")
  invisible(x)
}

# Prints a summary of the items that have been read, the table itself is still saved in the object.
#' @export
print.asari_project_results <- function(x, ...) {
  invisible(list(...))
  cat("asariR project results.\n")
  cat("Project directory: ", x$project_dir, "\n", sep = "")
  if (!is.null(x$preferred)) cat("Preferred features: ", nrow(x$preferred), "\n", sep = "")
  if (!is.null(x$full)) cat("Full features: ", nrow(x$full), "\n", sep = "")
  invisible(x)
}
