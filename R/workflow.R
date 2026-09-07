# Corresponds to Python asari/workflow.py: project scheduling, sample registration and mass track batch extraction.
# This module is only responsible for process orchestration; the actual mzML reading and mass track algorithm are located in chromatograms.R.

# Read list or environment parameters, compatible with Python mutable dict and R named list.
.workflow_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && name %in% names(object)) return(object[[name]])
  default
}

# Modify the parameters; the environment is updated in place, and the list returns an updated copy for the R caller to continue using.
.workflow_set <- function(object, name, value) {
  if (is.environment(object)) {
    assign(name, value, envir = object)
    return(object)
  }
  object[[name]] <- value
  object
}

# Read or write to the registry by sample ID, which always retains a Python 0-based value.
.workflow_registry_get <- function(registry, sample_id) {
  key <- as.character(sample_id)
  if (!is.null(names(registry)) && key %in% names(registry)) {
    return(registry[[key]])
  }
  registry[[as.integer(sample_id) + 1L]]
}

.workflow_registry_set <- function(registry, sample_id, value) {
  key <- as.character(sample_id)
  if (!is.null(names(registry)) && key %in% names(registry)) {
    registry[[key]] <- value
  } else {
    registry[[as.integer(sample_id) + 1L]] <- value
  }
  registry
}

# Get the month, day, hour, minute, and second string consistent with Python time.localtime()[1:6].
.workflow_time_parts <- function(now = Sys.time()) {
  local <- as.POSIXlt(now)
  as.character(c(
    local$mon + 1L, local$mday, local$hour, local$min,
    as.integer(local$sec)
  ))
}

# Find mass track extractor from parameters injection or chromatograms.R.
.workflow_mass_track_extractor <- function(parameters) {
  callback <- .workflow_get(parameters, "extract_mass_tracks")
  if (is.function(callback)) return(callback)
  if (exists("extract_mass_tracks", mode = "function")) return(extract_mass_tracks)
  stop("extract_mass_tracks is unavailable; load chromatograms.R first.")
}

# Call mass2chem's equivalent m/z difference anchor point finder.
.workflow_find_anchor_pairs <- function(list_mass_tracks, parameters) {
  callback <- .workflow_get(parameters, "find_mzdiff_pairs")
  if (is.function(callback)) {
    return(callback(
      list_mass_tracks,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm
    ))
  }
  if (exists(".constructors_find_mzdiff_pairs", mode = "function")) {
    return(.constructors_find_mzdiff_pairs(
      list_mass_tracks,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm
    ))
  }
  stop("m/z anchor-pair finder is unavailable; load constructors.R first.")
}

# Write the single-sample intermediate object as Python pickle to ensure that both samples.R and Python can read it.
.workflow_write_intermediate <- function(value, path) {
  if (!exists(".experiment_write_pickle", mode = "function")) {
    stop("Python pickle writer is unavailable; load experiment.R first.")
  }
  .experiment_write_pickle(value, path)
}

# Corresponds to workflow_setup: establish the output directory, extract EIC in batches and construct ext_Experiment.
workflow_setup <- function(list_input_files, parameters) {
  sample_registry <- register_samples(list_input_files)

  # auto mode selects memory or ondisk by project size, the threshold is exactly the same as in Python.
  if (identical(parameters$database_mode, "auto")) {
    database_mode <- if (
        length(list_input_files) <= parameters$project_sample_number_small) {
      "memory"
    } else {
      "ondisk"
    }
    parameters <- .workflow_set(parameters, "database_mode", database_mode)
  }

  time_parts <- .workflow_time_parts()
  parameters <- .workflow_set(
    parameters, "time_stamp", paste(time_parts, collapse = ":")
  )
  parameters <- create_export_folders(
    parameters, paste(time_parts, collapse = "")
  )
  shared_dict <- batch_EIC_from_samples_(sample_registry, parameters)
  if (length(shared_dict) == 0L) {
    stop("No data was processed, check the input files.", call. = FALSE)
  }

  for (sample_id in seq_along(sample_registry) - 1L) {
    sample <- .workflow_registry_get(sample_registry, sample_id)
    result <- .workflow_registry_get(shared_dict, sample_id)
    if (is.null(result)) {
      stop("Missing extraction result for sample ID ", sample_id, call. = FALSE)
    }

    # The current source code of Python misplaces the assignment of the last three tuple fields; here press single_sample_EICs_
    # The actual return order restores acquisition_time, sample_data, and sparsified to ensure that the process can be run.
    field_names <- c(
      "status:mzml_parsing", "status:eic", "data_location",
      "max_scan_number", "list_scan_numbers", "list_retention_time",
      "track_mzs", "number_anchor_mz_pairs", "anchor_mz_pairs",
      "acquisition_time", "sample_data", "sparsified"
    )
    for (ii in seq_along(field_names)) sample[[field_names[[ii]]]] <- result[[ii]]
    sample$name <- sub("\\.mzML", "", basename(sample$input_file), fixed = FALSE)
    sample_registry <- .workflow_registry_set(sample_registry, sample_id, sample)
  }

  ext_Experiment__init__(sample_registry, parameters)
}

# Corresponds to process_project: perform experimental processing, export and clean disk intermediate files according to workflow.
process_project <- function(list_input_files, parameters = default_parameters()) {
  experiment <- workflow_setup(list_input_files, parameters)
  workflow <- experiment$parameters$workflow
  message("Processing Experiment Using ", workflow, " Workflow...")

  if (workflow %in% c("LC", "GC")) {
    experiment$process_all_LC()
  } else if (identical(workflow, "DIMS")) {
    experiment$process_all_DIMS()
  } else if (identical(workflow, "LCMSMS")) {
    experiment$process_all_LCMSMS()
  } else {
    message("Workflow not supported. Error 100.")
  }
  experiment$export_all(anno = isTRUE(experiment$parameters$anno))

  if (!identical(experiment$parameters$database_mode, "memory")) {
    remove_intermediate_pickles(experiment$parameters)
  }
  invisible(NULL)
}

# Corresponds to read_project_dir: returns the file path whose name contains file_pattern in the directory.
read_project_dir <- function(directory, file_pattern = ".mzML") {
  message("Working on ~~ ", directory, " ~~ \n")
  entries <- list.files(directory, all.files = FALSE, no.. = TRUE)
  entries <- entries[grepl(file_pattern, entries, fixed = TRUE)]
  file.path(directory, entries)
}

# Corresponds to read_project_file: Read absolute path containing file_pattern from a text file.
read_project_file <- function(project_file, file_pattern = ".mzML") {
  message("Working on ~~ ", project_file, " ~~ \n")
  lines <- readLines(project_file, warn = FALSE)
  selected <- trimws(lines[grepl(file_pattern, lines, fixed = TRUE)])
  vapply(selected, function(path) {
    expanded <- path.expand(path)
    is_absolute <- grepl("^/", expanded) ||
      grepl("^[A-Za-z]:[/\\\\]", expanded)
    candidate <- if (is_absolute) expanded else file.path(getwd(), expanded)
    normalizePath(candidate, mustWork = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

# Corresponds to register_samples: Python 0-based sample_id is assigned to the input file in turn.
register_samples <- function(list_input_files) {
  registry <- lapply(seq_along(list_input_files), function(ii) {
    list(sample_id = ii - 1L, input_file = list_input_files[[ii]])
  })
  names(registry) <- as.character(seq_along(list_input_files) - 1L)
  registry
}

# Corresponds to create_export_folders: Create project, export and pickle directories and update parameters.
create_export_folders <- function(parameters, time_stamp = NULL) {
  if (grepl(parameters$project_name, parameters$outdir, fixed = TRUE)) {
    message("Export folders already exist, will not overwrite.")
    return(parameters)
  }

  suffix <- if (is.null(time_stamp)) {
    parameters$time_stamp_for_dir
  } else {
    time_stamp
  }
  outdir <- paste(parameters$outdir, parameters$project_name, suffix, sep = "_")
  parameters <- .workflow_set(parameters, "outdir", outdir)

  # Python os.makedirs throws an error when the target already exists; R also does not silently overwrite existing projects.
  if (!dir.create(outdir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Failed to create output directory: ", outdir, call. = FALSE)
  }
  export_dir <- file.path(outdir, "export")
  if (!dir.create(export_dir, showWarnings = FALSE) && !dir.exists(export_dir)) {
    stop("Failed to create export directory: ", export_dir, call. = FALSE)
  }
  pickle_dir <- file.path(outdir, "pickle")
  if (!dir.create(pickle_dir, showWarnings = FALSE) && !dir.exists(pickle_dir)) {
    stop("Failed to create pickle directory: ", pickle_dir, call. = FALSE)
  }
  parameters <- .workflow_set(parameters, "export_outdir", export_dir)
  parameters <- .workflow_set(parameters, "tmp_pickle_dir", pickle_dir)
  parameters
}

# Corresponds to remove_intermediate_pickles: delete temporary files and empty directories in the project pickle directory.
remove_intermediate_pickles <- function(parameters) {
  if (!is.null(parameters$reuse_intermediates)) {
    stop("Cannot remove when reuse_intermediates is set.", call. = FALSE)
  }
  message("Removing temporary pickle files...")
  pickle_dir <- parameters$tmp_pickle_dir
  entries <- list.files(pickle_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  for (entry in entries) {
    if (dir.exists(entry)) {
      stop("Intermediate directory contains a subdirectory: ", entry)
    }
    if (!file.remove(entry)) stop("Failed to remove intermediate: ", entry)
  }
  # The files have been confirmed and deleted one by one before, here only the exact pickle_dir that is empty is removed.
  removal_status <- unlink(pickle_dir, recursive = TRUE, force = FALSE)
  if (removal_status != 0L) {
    message("Failed to remove directory ", pickle_dir, ".")
  }
  invisible(NULL)
}

# Corresponds to make_iter_parameters: generate extraction task tuple for each sample.
make_iter_parameters <- function(sample_registry, parameters) {
  lapply(sample_registry, function(sample) {
    sample_name <- sub("\\.mzML", "", basename(sample$input_file))
    outfile <- file.path(
      parameters$outdir, "pickle", paste0(sample_name, ".pickle")
    )
    list(sample$sample_id, sample$input_file, outfile, parameters)
  })
}

# Corresponds to batch_EIC_from_samples_: Execute single sample mass track extraction in batches and merge the results.
batch_EIC_from_samples_ <- function(sample_registry, parameters) {
  jobs <- make_iter_parameters(sample_registry, parameters)
  workers <- parameters$multicores
  use_parallel <- .Platform$OS.type != "windows" && is.numeric(workers) &&
    length(workers) == 1L && workers > 1L
  results <- if (use_parallel) {
    parallel::mclapply(jobs, single_sample_EICs_, mc.cores = as.integer(workers))
  } else {
    lapply(jobs, single_sample_EICs_)
  }
  merged <- list()
  for (result in results) {
    for (name in names(result)) merged[[name]] <- result[[name]]
  }
  merged
}

# Corresponds to single_sample_EICs_: extract mass track, anchor point and scan information of an mzML.
single_sample_EICs_ <- function(job) {
  sample_id <- job[[1L]]
  infile <- job[[2L]]
  outfile <- job[[3L]]
  parameters <- job[[4L]]

  tryCatch({
    extractor <- .workflow_mass_track_extractor(parameters)
    xdict <- extractor(
      infile,
      mz_tolerance_ppm = parameters$mz_tolerance_ppm,
      min_intensity = parameters$min_intensity_threshold,
      min_timepoints = parameters$min_timepoints,
      min_peak_height = parameters$min_peak_height
    )
    tracks <- .workflow_get(xdict, "tracks", list())
    list_mass_tracks <- lapply(seq_along(tracks), function(ii) {
      list(id_number = ii - 1L, mz = tracks[[ii]][[1L]],
           intensity = tracks[[ii]][[2L]])
    })
    message(
      "Extracted ", basename(infile), " to ", length(tracks), " mass tracks."
    )
    anchor_pairs <- .workflow_find_anchor_pairs(list_mass_tracks, parameters)
    rt_numbers <- .workflow_get(xdict, "rt_numbers")
    if (is.null(rt_numbers) || length(rt_numbers) == 0L) {
      stop("No MS1 scan number was extracted.")
    }
    new <- list(
      sample_id = sample_id,
      input_file = infile,
      ion_mode = parameters$mode,
      list_mass_tracks = list_mass_tracks,
      anchor_mz_pairs = anchor_pairs,
      number_anchor_mz_pairs = length(anchor_pairs),
      track_mzs = lapply(seq_along(tracks), function(ii) {
        c(tracks[[ii]][[1L]], ii - 1L)
      }),
      ms2_spectra = .workflow_get(xdict, "ms2_spectra", list()),
      max_scan_number = max(rt_numbers),
      acquisition_time = .workflow_get(xdict, "acquisition_time")
    )

    data_filepath <- outfile
    if (identical(parameters$database_mode, "ondisk")) {
      .workflow_write_intermediate(new, outfile)
    }
    tuple <- list(
      "passed", "passed", data_filepath, new$max_scan_number,
      rt_numbers, .workflow_get(xdict, "rt_times"), new$track_mzs,
      new$number_anchor_mz_pairs, anchor_pairs, new$acquisition_time,
      if (identical(parameters$database_mode, "memory")) new else list(),
      parameters$compress
    )
    result <- list(tuple)
    names(result) <- as.character(sample_id)
    result
  }, error = function(error) {
    message("Failed to extract ", basename(infile), ": ", conditionMessage(error))
    tuple <- list(
      "failed", "failed", NULL, NULL, NULL, NULL, NULL, NULL,
      NULL, NULL, NULL, parameters$compress
    )
    result <- list(tuple)
    names(result) <- as.character(sample_id)
    result
  })
}

# Corresponds to process_xics: only extract and save XIC, do not run Experiment/FeatureTable.
process_xics <- function(list_input_files, parameters) {
  parameters <- .workflow_set(parameters, "database_mode", "ondisk")
  time_stamp <- paste(.workflow_time_parts(), collapse = "")
  parameters <- create_export_folders(parameters, time_stamp)

  # The current source code of Python is mistakenly written as register_samples(register_samples(files));
  # R uses one-time registration to avoid treating sample IDs as mzML paths.
  batch_EIC_from_samples_(register_samples(list_input_files), parameters)
  message(
    "XICs were stored as pickle objects under ",
    file.path(parameters$outdir, "pickle")
  )
  invisible(NULL)
}

# Corresponds to get_mz_list: skip the header and read the target m/z separated by tab/comma in the first column.
get_mz_list <- function(infile) {
  lines <- readLines(infile, warn = FALSE)
  if (length(lines) <= 1L) return(numeric())
  values <- vapply(lines[-1L], function(line) {
    first_tab_field <- strsplit(line, "\t", fixed = TRUE)[[1L]][[1L]]
    first_comma_field <- strsplit(first_tab_field, ",", fixed = TRUE)[[1L]][[1L]]
    value <- suppressWarnings(as.numeric(first_comma_field))
    if (is.na(value)) stop("Invalid m/z value: ", first_comma_field)
    value
  }, numeric(1))
  unname(values)
}
