# R package entry point: validate explicit inputs, organize parameters, and call the core asari workflow.

# Check out the R dependencies necessary to actually process mzML and write out project records.
.asari_check_runtime_dependencies <- function() {
  required <- c("mzR", "jsonlite")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      ". Install mzR with BiocManager and jsonlite with install.packages().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Collect mzML from files or directories explicitly given by the user; hidden default data locations are not searched.
.asari_collect_input_files <- function(input, recursive = FALSE) {
  if (!is.character(input) || length(input) == 0L || anyNA(input) ||
      any(!nzchar(trimws(input)))) {
    stop("input must contain at least one explicit mzML file or directory.", call. = FALSE)
  }
  if (!is.logical(recursive) || length(recursive) != 1L || is.na(recursive)) {
    stop("recursive must be TRUE or FALSE.", call. = FALSE)
  }

  # Expand the path item by item, the directory only reads standard.mzML files, and the extension of a single file is strictly checked.
  collected <- character()
  for (item in path.expand(input)) {
    if (dir.exists(item)) {
      collected <- c(
        collected,
        list.files(
          item,
          pattern = "[.]mzML$",
          full.names = TRUE,
          recursive = recursive,
          ignore.case = FALSE
        )
      )
    } else if (file.exists(item)) {
      if (!grepl("[.]mzML$", basename(item))) {
        stop("Input file is not a .mzML file: ", item, call. = FALSE)
      }
      collected <- c(collected, item)
    } else {
      stop("Input path does not exist: ", item, call. = FALSE)
    }
  }

  # Normalize and sort to ensure that the order of samples from the same input is stable across different runs.
  collected <- sort(unique(normalizePath(collected, mustWork = TRUE)))
  if (length(collected) == 0L) {
    stop("No .mzML files were found in the supplied input path(s).", call. = FALSE)
  }

  # Asari uses file names to register samples. Duplicate file names will make the output columns difficult to distinguish, so an error is reported in advance.
  duplicated_names <- unique(basename(collected)[duplicated(basename(collected))])
  if (length(duplicated_names) > 0L) {
    stop(
      "Duplicate mzML file name(s) are not supported: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }
  collected
}

# Merge user overrides into a new default parameter to avoid modifying module-level PARAMETERS objects.
.asari_merge_parameters <- function(overrides = list()) {
  if (!is.list(overrides)) {
    stop("parameters must be a named list.", call. = FALSE)
  }
  if (length(overrides) > 0L &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("Every parameters entry must have a name.", call. = FALSE)
  }
  parameters <- default_parameters()
  if (length(overrides) > 0L) parameters[names(overrides)] <- overrides
  parameters
}

# Verify the most commonly used parameters in the public entry and write back the internal parameter list.
.asari_prepare_parameters <- function(
    parameters, output, project_name, mode = NULL, ppm = NULL,
    multicores = NULL, rt_align = NULL, database_mode = NULL) {
  if (!is.character(project_name) || length(project_name) != 1L ||
      is.na(project_name) || !nzchar(trimws(project_name)) ||
      grepl("[/\\\\]", project_name)) {
    stop("project_name must be one non-empty name without path separators.", call. = FALSE)
  }
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(trimws(output))) {
    stop("output must be one explicit directory path.", call. = FALSE)
  }

  # Create a user-specified parent directory; each run will still generate an independent timestamp project directory by the workflow.
  output <- path.expand(output)
  if (file.exists(output) && !dir.exists(output)) {
    stop("output exists but is not a directory: ", output, call. = FALSE)
  }
  if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create output directory: ", output, call. = FALSE)
  }
  output <- normalizePath(output, mustWork = TRUE)

  # Explicit parameters take precedence over the parameters list; the user list or asari defaults are retained when not given explicitly.
  if (!is.null(mode)) parameters$mode <- mode
  if (!parameters$mode %in% c("pos", "neg")) {
    stop("mode must be either 'pos' or 'neg'.", call. = FALSE)
  }
  if (!is.null(ppm)) parameters$mz_tolerance_ppm <- ppm
  if (!is.numeric(parameters$mz_tolerance_ppm) ||
      length(parameters$mz_tolerance_ppm) != 1L ||
      is.na(parameters$mz_tolerance_ppm) || parameters$mz_tolerance_ppm <= 0) {
    stop("ppm must be one positive number.", call. = FALSE)
  }
  if (!is.null(multicores)) parameters$multicores <- multicores
  if (!is.numeric(parameters$multicores) || length(parameters$multicores) != 1L ||
      is.na(parameters$multicores) || parameters$multicores < 1 ||
      parameters$multicores != as.integer(parameters$multicores)) {
    stop("multicores must be one positive integer.", call. = FALSE)
  }
  parameters$multicores <- as.integer(parameters$multicores)
  if (!is.null(rt_align)) parameters$rt_align_on <- rt_align
  if (!is.logical(parameters$rt_align_on) || length(parameters$rt_align_on) != 1L ||
      is.na(parameters$rt_align_on)) {
    stop("rt_align must be TRUE or FALSE.", call. = FALSE)
  }

  # The new convenience entry uses memory mode by default, avoiding the need for ordinary LC processing to rely on Python pickle conversion.
  if (is.null(database_mode)) {
    database_mode <- if ("database_mode" %in% names(parameters)) {
      parameters$database_mode
    } else {
      "memory"
    }
  }
  if (!database_mode %in% c("memory", "ondisk", "auto")) {
    stop("database_mode must be 'memory', 'ondisk', or 'auto'.", call. = FALSE)
  }
  parameters$database_mode <- database_mode
  parameters$project_name <- project_name
  parameters$outdir <- file.path(output, "run")
  parameters
}

# Compute the derived peak-detection thresholds in the order of Python command line entry, while retaining the advanced coverage value explicitly given by the user.
.asari_finalize_peak_parameters <- function(parameters, input_files, override_names = character()) {
  if (!is.logical(parameters$autoheight) || length(parameters$autoheight) != 1L ||
      is.na(parameters$autoheight)) {
    stop("parameters$autoheight must be TRUE or FALSE.", call. = FALSE)
  }

  # Automatic peak height directly uses the explicit files that have already been collected, avoiding any default directories to be scanned again.
  if (isTRUE(parameters$autoheight)) {
    if (!exists("estimate_min_peak_height", mode = "function", inherits = TRUE)) {
      stop(
        "autoheight requires the optional analyze module; set autoheight = FALSE ",
        "or install the complete asariR source tree.",
        call. = FALSE
      )
    }
    estimator <- get("estimate_min_peak_height", mode = "function", inherits = TRUE)
    parameters$min_peak_height <- estimator(input_files, parameters)
    parameters$min_intensity_threshold <- parameters$min_peak_height / 10
  }

  # The Python entry calculates both values based on the final peak height; the explicit value in the high-level list has highest priority.
  if (!("min_prominence_threshold" %in% override_names)) {
    parameters$min_prominence_threshold <- as.integer(0.33 * parameters$min_peak_height)
  }
  if (!("cal_min_peak_height" %in% override_names)) {
    parameters$cal_min_peak_height <- 10 * parameters$min_peak_height
  }

  # Check all critical thresholds at once before entering time-consuming processing, and give users clear parameter names.
  positive_fields <- c(
    "min_peak_height", "min_prominence_threshold",
    "cal_min_peak_height", "min_intensity_threshold"
  )
  for (name in positive_fields) {
    value <- parameters[[name]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
      stop("parameters$", name, " must be one positive number.", call. = FALSE)
    }
  }
  parameters
}

#' Use asariR to process centroided mzML data
#'
#' This is the main entry point for users. The function receives an mzML file or directory explicitly specified by the user,
#' Runs the translated asari process and returns the location of the main results table.
#'
#' @param input One or more explicit `.mzML` files, or a directory containing `.mzML` files.
#' @param output The parent directory that holds the timestamped project directory.
#' @param project_name Short project name for the output directory name.
#' @param mode Optional ion mode can only be `"pos"` or `"neg"`.
#' @param ppm Optional positive m/z tolerance in ppm.
#' @param multicores Optional positive integer number of worker processes.
#' @param rt_align Optional logic value, controls whether to perform retention time alignment.
#' @param database_mode Data storage mode. The convenience entry point uses `"memory"` by default;
#' `"ondisk"` and some `"auto"` require a Python interpreter to be compatible with pickle to run.
#' @param recursive Whether to search input directories recursively.
#' @param parameters Advanced parameter override list from [asari_default_parameters()];
#' Explicit parameters have higher precedence.
#'
#' @return An `asari_result` object containing the project directory, feature-table paths, input files, and validated parameters.
#' @export
asari_process <- function(
    input,
    output = "asari-output",
    project_name = "asari_project",
    mode = NULL,
    ppm = NULL,
    multicores = NULL,
    rt_align = NULL,
    database_mode = NULL,
    recursive = FALSE,
    parameters = list()) {
  # Complete lightweight input validation first, and then check the software dependencies required to process real mzML.
  input_files <- .asari_collect_input_files(input, recursive = recursive)
  .asari_check_runtime_dependencies()

  # Incorporate advanced parameters, then apply the common explicit parameters exposed by the public entry point.
  effective <- .asari_merge_parameters(parameters)
  if (is.null(database_mode) && !("database_mode" %in% names(parameters))) {
    database_mode <- "memory"
  }
  effective <- .asari_prepare_parameters(
    effective,
    output = output,
    project_name = project_name,
    mode = mode,
    ppm = ppm,
    multicores = multicores,
    rt_align = rt_align,
    database_mode = database_mode
  )
  effective <- .asari_finalize_peak_parameters(
    effective,
    input_files = input_files,
    override_names = names(parameters)
  )

  # The disk intermediate mode requires the Python standard library to write a pickle file compatible with the original asari and give clear errors in advance.
  if (!identical(effective$database_mode, "memory")) {
    .samples_find_python()
  }

  # Record the directories present before the run, and the difference set is used to accurately locate the timestamp item generated this time.
  project_pattern <- paste0(
    effective$outdir, "_", effective$project_name, "_*"
  )
  existing_projects <- Sys.glob(project_pattern)
  process_project(input_files, effective)
  generated_projects <- setdiff(Sys.glob(project_pattern), existing_projects)
  if (length(generated_projects) != 1L) {
    stop(
      "Processing finished but the new project directory could not be identified.",
      call. = FALSE
    )
  }

  # Return the most commonly used result paths, while retaining input and valid parameters for review and reproduction.
  project_dir <- normalizePath(generated_projects[[1L]], mustWork = TRUE)
  # The output file name allows advanced parameter overrides, so the return path must use the final effective name.
  feature_table_name <- effective$output_feature_table
  result <- list(
    project_dir = project_dir,
    preferred_feature_table = file.path(
      project_dir, paste0("preferred_", feature_table_name)
    ),
    full_feature_table = file.path(
      project_dir, "export", paste0("full_", feature_table_name)
    ),
    input_files = input_files,
    parameters = effective
  )
  class(result) <- "asari_result"
  print(result)
  invisible(result)
}

#' Returns a separate copy of asariR's default parameters
#'
#' @return A named list that can be safely modified.
#' @export
asari_default_parameters <- function() {
  default_parameters()
}

# Print a concise summary of results to avoid users looking for the timestamp output directory themselves.
#' @export
print.asari_result <- function(x, ...) {
  # Accepts and ignores extra print parameters to conform to R's S3 print method convention.
  dots <- list(...)
  invisible(dots)
  cat("asariR processing completed.\n")
  cat("Project directory: ", x$project_dir, "\n", sep = "")
  cat("Preferred feature table: ", x$preferred_feature_table, "\n", sep = "")
  cat("Full feature table: ", x$full_feature_table, "\n", sep = "")
  invisible(x)
}
