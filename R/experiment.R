# Corresponds to Python asari/experiment.py: experimental-level sample management, processing, annotation and result export.
# Python's mutable ext_Experiment class is implemented in R using a class-tagged environment.

# Read a list or environment field; this helper only handles R/Python object structure differences.
.experiment_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && name %in% names(object)) return(object[[name]])
  default
}

# Read the registry by sample ID, compatible with named IDs and Python 0-based numeric IDs.
.experiment_registry_get <- function(registry, sample_id) {
  key <- as.character(sample_id)
  if (is.environment(registry) && exists(key, registry, inherits = FALSE)) {
    return(get(key, registry, inherits = FALSE))
  }
  if (!is.list(registry)) stop("sample_registry must be a list or environment.")
  if (!is.null(names(registry)) && key %in% names(registry)) {
    return(registry[[key]])
  }
  if (is.numeric(sample_id)) return(registry[[as.integer(sample_id) + 1L]])
  stop("Unknown sample id: ", key, call. = FALSE)
}

# Returns the Python dict keys of the registry; unnamed lists use 0-based position.
.experiment_registry_ids <- function(registry) {
  if (is.environment(registry)) return(ls(registry, all.names = TRUE))
  if (length(registry) == 0L) return(vector(mode = "list", length = 0L))
  if (!is.null(names(registry)) && all(nzchar(names(registry)))) {
    registry_names <- names(registry)
    # The Python key of asari sample_registry is an integer; the R list name will automatically convert it to a character.
    if (all(grepl("^-?[0-9]+$", registry_names))) {
      return(as.list(as.integer(registry_names)))
    }
    return(as.list(registry_names))
  }
  as.list(seq_along(registry) - 1L)
}

# Calls a Python-style method on environment/list and preserves mutable side effects on the object.
.experiment_call <- function(object, method, ...) {
  callback <- .experiment_get(object, method)
  if (!is.function(callback)) {
    stop("Object does not provide method: ", method, call. = FALSE)
  }
  callback(...)
}

# Modify a value in parameters and rebind to the experimental object.
.experiment_set_parameter <- function(self, name, value) {
  parameters <- self$parameters
  parameters[[name]] <- value
  self$parameters <- parameters
  invisible(value)
}

# Write JSON, the function corresponds to Python json.dump(..., cls=NpEncoder).
.experiment_write_json <- function(value, path, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("JSON export requires the optional jsonlite package.", call. = FALSE)
  }
  text <- jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    pretty = pretty,
    digits = NA,
    dataframe = "rows"
  )
  writeLines(enc2utf8(text), path, useBytes = TRUE)
  invisible(path)
}

# Write a real Python pickle so that the asari Python tool continues to read R-exported CMAPs.
.experiment_write_pickle <- function(value, path) {
  if (!exists(".samples_find_python", mode = "function")) {
    stop("Python pickle export requires samples.R to be loaded.", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Python pickle export requires the optional jsonlite package.")
  }
  python <- .samples_find_python()
  json_path <- tempfile(fileext = ".json")
  script_path <- tempfile(fileext = ".py")
  on.exit(unlink(c(json_path, script_path)), add = TRUE)
  .experiment_write_json(value, json_path, pretty = FALSE)
  writeLines(c(
    "import json, pickle, sys",
    "source, output = sys.argv[1:3]",
    "with open(source, 'r', encoding='utf-8') as handle: value = json.load(handle)",
    "with open(output, 'wb') as handle: pickle.dump(value, handle, pickle.HIGHEST_PROTOCOL)"
  ), script_path)
  output <- system2(
    python,
    shQuote(c(script_path, json_path, path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Failed to write Python pickle: ", paste(output, collapse = "\n"))
  }
  invisible(path)
}

# experiment.py imports the annotation pattern and README text from default_parameters.py.
.experiment_adduct_patterns_pos <- list(
  c(21.98194, "Na/H"), c(41.026549, "ACN"),
  c(67.987424, "NaCOOH"), c(37.955882, "K/H")
)
.experiment_adduct_patterns_neg <- list(
  c(21.98194, "Na/H"), c(67.987424, "NaCOOH"),
  c(82.0030, "C2HF3"), c(1.99566, "F <-> OH"),
  c(41.026549, "ACN"), c(37.955882, "K/H")
)
.experiment_isotope_patterns <- list(
  list(1.003355, "13C/12C", c(0, 0.8)),
  list(2.00671, "13C/12C*2", c(0, 0.8)),
  list(3.010065, "13C/12C*3", c(0, 0.8)),
  list(1.9970, "37Cl/35Cl", c(0.1, 0.8))
)
.experiment_extended_adducts <- list(
  c(1.0078, "H"), c(-1.0078, "-H"),
  c(10.991, "Na/H, double charged"), c(0.5017, "13C/12C, double charged"),
  c(117.02655, "-NH3"), c(17.02655, "NH3"),
  c(-18.0106, "-H2O"), c(18.0106, "H2O"), c(18.033823, "NH4"),
  c(27.01089904, "HCN"), c(37.94694, "Ca/H2"), c(32.026215, "MeOH"),
  c(43.96389, "Na2/H2"), c(67.987424, "NaCOOH"),
  c(83.961361, "KCOOH"), c(97.96737927, "H2SO4"),
  c(97.97689507, "H3PO4")
)
.experiment_readme <- paste0(
  "\nThe recommended feature table is `preferred_Feature_table.tsv`. \n\n",
  "All peaks are kept in `export/full_Feature_table.tsv` \n",
  "if they meet signal (snr) and shape standards \n",
  "(part of input parameters but default values are fine for most people). \n",
  "The filtering decisions are left to end users.\n\n",
  "Annotation is in JSON (`Annotated_empricalCompounds.json`) \n",
  "and in tab delimited text (`Feature_annotation.tsv`).\n\n",
  "The processing parameters and history are in `project.json`.\n\n",
  "Please refer to https://github.com/shuzhao-li-lab/asari for details, \n",
  "report bugs or request features.\n"
)

# Corresponds to ext_Experiment.__init__: Creates experimental mutable class objects.
ext_Experiment__init__ <- function(sample_registry, parameters) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("ext_Experiment", "environment")
  self$sample_registry <- sample_registry
  self$parameters <- parameters

  # Python's all_samples and all_sample_instances point to the same list; active binding retains aliases.
  self$`.all_samples_store` <- list()
  shared_samples <- function(value) {
    if (!missing(value)) self$`.all_samples_store` <- value
    self$`.all_samples_store`
  }
  makeActiveBinding("all_samples", shared_samples, self)
  makeActiveBinding("all_sample_instances", shared_samples, self)

  # First install the methods required in the construction phase, and then calculate the number of samples, number of scans and reference sample.
  self$get_valid_sample_ids <- function() ext_Experiment_get_valid_sample_ids(self)
  self$get_max_scan_number <- function(sample_registry) {
    ext_Experiment_get_max_scan_number(self, sample_registry)
  }
  self$get_reference_sample_id <- function() {
    ext_Experiment_get_reference_sample_id(self)
  }
  self$valid_sample_ids <- ext_Experiment_get_valid_sample_ids(self)
  self$number_of_samples <- length(self$valid_sample_ids)
  self$number_scans <- ext_Experiment_get_max_scan_number(self, sample_registry)

  self$output_dir <- parameters$outdir
  self$mode <- parameters$mode
  self$mz_tolerance_ppm <- parameters$mz_tolerance_ppm
  self$check_isotope_ratio <- parameters$check_isotope_ratio
  self$database_mode <- parameters$database_mode
  self$reference_sample_id <- ext_Experiment_get_reference_sample_id(self)

  # constructors.py imports SimpleSample directly in Python; R maintains the same dependency via factory.
  if (exists("SimpleSample__init__", mode = "function")) {
    self$sample_factory <- function(...) {
      sample <- SimpleSample__init__(...)
      # SimpleSample holds Python 0-based track IDs; R index bounds for constructors use 1-based.
      landmarks <- sample$`_mz_landmarks_`
      if (!is.null(landmarks) && length(landmarks) > 0L) {
        sample$`_mz_landmarks_` <- as.integer(landmarks) + 1L
      }
      sample
    }
  }

  # Install the remaining Python instance methods for the environment.
  self$determine_acquisition_order <- function() {
    ext_Experiment_determine_acquisition_order(self)
  }
  self$process_all_LC <- function() ext_Experiment_process_all_LC(self)
  self$process_all_DIMS <- function() ext_Experiment_process_all_DIMS(self)
  self$process_all_LCMSMS <- function() ext_Experiment_process_all_LCMSMS(self)
  self$export_all <- function(anno = FALSE) ext_Experiment_export_all(self, anno)
  self$annotate <- function() ext_Experiment_annotate(self)
  self$generate_qc_plot_pdf <- function(outfile = "qc_plot.pdf") {
    ext_Experiment_generate_qc_plot_pdf(self, outfile)
  }
  self$export_CMAP_pickle <- function() ext_Experiment_export_CMAP_pickle(self)
  self$load_annotation_db <- function(src = "hmdb4") {
    ext_Experiment_load_annotation_db(self, src)
  }
  self$db_mass_calibrate <- function(required_calibrate_threshold = 0.000002) {
    ext_Experiment_db_mass_calibrate(self, required_calibrate_threshold)
  }
  self$append_orphans_to_epmCpds <- function(dict_empCpds) {
    ext_Experiment_append_orphans_to_epmCpds(self, dict_empCpds)
  }
  self$export_peak_annotation <- function(dict_empCpds, KCD, prefix) {
    ext_Experiment_export_peak_annotation(self, dict_empCpds, KCD, prefix)
  }
  self$select_unique_compound_features <- function(dict_empCpds) {
    ext_Experiment_select_unique_compound_features(self, dict_empCpds)
  }
  self$export_feature_tables <- function(
      `_snr` = 2, `_peak_shape` = 0.7, `_cSelectivity` = 0.7) {
    ext_Experiment_export_feature_tables(
      self, `_snr`, `_peak_shape`, `_cSelectivity`
    )
  }
  self$export_log <- function() ext_Experiment_export_log(self)
  self$export_readme <- function() ext_Experiment_export_readme(self)
  self
}

# Corresponds to ext_Experiment.get_reference_sample_id: select reference sample according to the specified file or anchor number.
ext_Experiment_get_reference_sample_id <- function(self) {
  reference <- self$parameters$reference
  ids <- .experiment_registry_ids(self$sample_registry)
  if (!is.null(reference) && length(reference) > 0L && nzchar(reference)) {
    for (sample_id in ids) {
      registry <- .experiment_registry_get(self$sample_registry, sample_id)
      input_file <- .experiment_get(registry, "input_file")
      if (basename(reference) == basename(input_file) ||
          paste0(basename(reference), ".mzML") == basename(input_file)) {
        return(sample_id)
      }
    }
    return(NULL)
  }

  if (length(ids) == 0L) return(NULL)
  candidates <- lapply(ids, function(sample_id) {
    registry <- .experiment_registry_get(self$sample_registry, sample_id)
    anchors <- .experiment_get(registry, "number_anchor_mz_pairs")
    if (is.null(anchors)) return(NULL)
    list(
      anchors = anchors,
      sample_id = .experiment_get(registry, "sample_id")
    )
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (length(candidates) == 0L) stop("No sample has number_anchor_mz_pairs.")
  ordering <- order(
    vapply(candidates, `[[`, numeric(1), "anchors"),
    vapply(candidates, function(x) as.numeric(x$sample_id), numeric(1)),
    decreasing = TRUE
  )
  reference_id <- candidates[[ordering[[1L]]]]$sample_id
  reference_registry <- .experiment_registry_get(
    self$sample_registry, reference_id
  )
  .experiment_set_parameter(
    self, "reference", .experiment_get(reference_registry, "input_file")
  )
  message(
    "\n    The reference sample is:\n    ||* ",
    .experiment_get(reference_registry, "name"), " *||\n"
  )
  message(
    "Max reference retention time is ",
    sprintf("%4.2f", max(.experiment_get(reference_registry, "list_retention_time"))),
    " at scan number ", .experiment_get(reference_registry, "max_scan_number"),
    ".\n"
  )
  reference_id
}

# Corresponds to ext_Experiment.get_valid_sample_ids: only retain samples with successful EIC extraction.
ext_Experiment_get_valid_sample_ids <- function(self) {
  ids <- .experiment_registry_ids(self$sample_registry)
  Filter(function(sample_id) {
    identical(
      .experiment_get(
        .experiment_registry_get(self$sample_registry, sample_id), "status:eic"
      ),
      "passed"
    )
  }, ids)
}

# Corresponds to ext_Experiment.determine_acquisition_order: priority is based on acquisition time, and in case of failure, based on ID.
ext_Experiment_determine_acquisition_order <- function(self) {
  ids <- .experiment_registry_ids(self$sample_registry)
  pairs <- lapply(ids, function(sample_id) {
    list(
      sample_id,
      .experiment_get(
        .experiment_registry_get(self$sample_registry, sample_id),
        "acquisition_time"
      )
    )
  })
  tryCatch({
    times <- lapply(pairs, `[[`, 2L)
    if (any(vapply(times, is.null, logical(1)))) stop("incomparable time")
    pairs[order(vapply(times, as.character, character(1)))]
  }, error = function(error) {
    numeric_ids <- suppressWarnings(as.numeric(vapply(
      pairs, function(pair) as.character(pair[[1L]]), character(1)
    )))
    if (all(!is.na(numeric_ids))) pairs[order(numeric_ids)] else {
      pairs[order(vapply(pairs, function(pair) {
        as.character(pair[[1L]])
      }, character(1)))]
    }
  })
}

# Corresponds to ext_Experiment.get_max_scan_number: the maximum 0-based scan plus one to get the total length.
ext_Experiment_get_max_scan_number <- function(self, sample_registry) {
  if (length(.experiment_registry_ids(sample_registry)) == 0L) return(NULL)
  if (length(self$valid_sample_ids) == 0L) stop("No valid sample is available.")
  max(vapply(self$valid_sample_ids, function(sample_id) {
    .experiment_get(
      .experiment_registry_get(sample_registry, sample_id), "max_scan_number"
    )
  }, numeric(1))) + 1L
}

# Create a CompositeMap; tests can provide a factory and call constructors.R during normal operation.
.experiment_make_composite_map <- function(self) {
  factory <- .experiment_get(self, "CompositeMap_factory")
  if (is.function(factory)) return(factory(self))
  CompositeMap__init__(self)
}

# Methods provided by the test or extension object are called first, otherwise the corresponding function of constructors.R is called.
.experiment_cmap_call <- function(cmap, method, ...) {
  callback <- .experiment_get(cmap, method)
  if (is.function(callback)) return(callback(...))
  function_name <- paste0("CompositeMap_", method)
  if (!exists(function_name, mode = "function")) {
    stop("Missing CompositeMap method: ", method, call. = FALSE)
  }
  get(function_name, mode = "function")(cmap, ...)
}

# Corresponds to ext_Experiment.process_all_LC: execute LC/GC MassGrid, alignment, synthesis and peak detection.
ext_Experiment_process_all_LC <- function(self) {
  self$CMAP <- .experiment_make_composite_map(self)
  .experiment_cmap_call(self$CMAP, "construct_mass_grid")
  if (!isTRUE(self$parameters$rt_align_on)) {
    .experiment_cmap_call(self$CMAP, "mock_rentention_alignment")
  }
  .experiment_cmap_call(self$CMAP, "build_composite_tracks")
  .experiment_cmap_call(self$CMAP, "global_peak_detection")
  invisible(NULL)
}

# Corresponds to ext_Experiment.process_all_DIMS: Only create MassGrid and extract the maximum intensity table.
ext_Experiment_process_all_DIMS <- function(self) {
  self$CMAP <- .experiment_make_composite_map(self)
  .experiment_cmap_call(self$CMAP, "construct_mass_grid")
  .experiment_cmap_call(self$CMAP, "get_DIMS_feature_table")
  invisible(NULL)
}

# Corresponds to ext_Experiment.process_all_LCMSMS: Python is currently also an empty placeholder.
ext_Experiment_process_all_LCMSMS <- function(self) {
  invisible(NULL)
}

# Corresponds to ext_Experiment.export_all: export the MassGrid, features, logs, and README for each workflow.
ext_Experiment_export_all <- function(self, anno = FALSE) {
  workflow <- self$parameters$workflow
  export_dir <- file.path(self$parameters$outdir, "export")
  if (workflow %in% c("LC", "GC")) {
    utils::write.csv(
      self$CMAP$MassGrid,
      file.path(export_dir, self$parameters$mass_grid_mapping),
      row.names = TRUE
    )
    if (isTRUE(anno)) {
      for (ii in seq_along(self$CMAP$FeatureList)) {
        self$CMAP$FeatureList[[ii]]$id <- as.character(
          self$CMAP$FeatureList[[ii]]$id_number
        )
      }
      ext_Experiment_export_CMAP_pickle(self)
      ext_Experiment_annotate(self)
    }
    ext_Experiment_export_feature_tables(self)
    ext_Experiment_export_log(self)
    ext_Experiment_export_readme(self)
  } else if (identical(workflow, "DIMS")) {
    utils::write.csv(
      self$CMAP$MassGrid,
      file.path(export_dir, self$parameters$mass_grid_mapping),
      row.names = TRUE
    )
    ext_Experiment_export_feature_tables(self)
    ext_Experiment_export_log(self)
  }
  invisible(NULL)
}

# Corresponds to ext_Experiment.annotate: completes three-stage annotation and file export through JMS compatible objects.
ext_Experiment_annotate <- function(self) {
  ext_Experiment_load_annotation_db(self)
  ext_Experiment_db_mass_calibrate(self)
  factory <- .experiment_get(self, "ExperimentalEcpdDatabase_factory")
  if (!is.function(factory)) {
    stop(
      "LC annotation requires an ExperimentalEcpdDatabase_factory adapter.",
      call. = FALSE
    )
  }
  EED <- factory(mode = self$mode, mz_tolerance_ppm = self$mz_tolerance_ppm,
                 rt_tolerance = 2)
  EED$adduct_patterns <- if (identical(self$mode, "pos")) {
    .experiment_adduct_patterns_pos
  } else {
    .experiment_adduct_patterns_neg
  }
  EED$isotope_search_patterns <- .experiment_isotope_patterns
  EED$extended_adducts <- .experiment_extended_adducts
  .experiment_call(EED, "build_from_list_peaks", self$CMAP$FeatureList)
  .experiment_call(EED, "extend_empCpd_annotation", self$KCD)
  .experiment_call(EED, "annotate_singletons", self$KCD)
  EED$dict_empCpds <- ext_Experiment_append_orphans_to_epmCpds(
    self, EED$dict_empCpds
  )
  ext_Experiment_export_peak_annotation(
    self, EED$dict_empCpds, self$KCD, "Feature_annotation"
  )
  if (length(.experiment_registry_ids(self$sample_registry)) > 0L) {
    ext_Experiment_select_unique_compound_features(self, EED$dict_empCpds)
  }
  .experiment_write_json(
    EED$dict_empCpds,
    file.path(self$parameters$outdir, "Annotated_empiricalCompounds.json")
  )
  .experiment_write_pickle(
    EED$dict_empCpds,
    file.path(self$parameters$outdir, "export", "epd.pickle")
  )
  invisible(NULL)
}

# Corresponds to ext_Experiment.generate_qc_plot_pdf: Call the QC plotter when available, otherwise skip.
ext_Experiment_generate_qc_plot_pdf <- function(self, outfile = "qc_plot.pdf") {
  plotter <- .experiment_get(self, "asari_qc_plot")
  if (!is.function(plotter)) {
    message("[QC plot] cannot import matplotlib, skipping.")
    return(invisible(NULL))
  }
  plotter(
    self$CMAP$FeatureTable,
    outfile = file.path(self$parameters$outdir, "export", outfile)
  )
  invisible(NULL)
}

# Corresponds to ext_Experiment.export_CMAP_pickle: Exports CMAP data that can be read visually by Python.
ext_Experiment_export_CMAP_pickle <- function(self) {
  rt_records <- lapply(self$all_samples, function(sample) {
    method <- .experiment_get(sample, "get_rt_calibration_records")
    if (is.function(method)) method() else {
      SimpleSample_get_rt_calibration_records(sample)
    }
  })
  export <- list(
    `_number_of_samples_` = self$CMAP$`_number_of_samples_`,
    rt_length = self$CMAP$rt_length,
    # Using a list preserves Python's JSON/pickle structure of a list even if there is only one element.
    rt_reference_landmarks = lapply(
      self$CMAP$good_reference_landmark_peaks, `[[`, "apex"
    ),
    rt_records = rt_records,
    dict_scan_rtime = as.list(self$CMAP$dict_scan_rtime),
    list_mass_tracks = self$CMAP$composite_mass_tracks,
    MassGrid = lapply(self$CMAP$MassGrid, function(column) {
      as.list(unname(column))
    })
  )
  .experiment_write_pickle(
    export,
    file.path(self$parameters$outdir, "export", "cmap.pickle")
  )
  invisible(NULL)
}

# Corresponds to ext_Experiment.load_annotation_db: Loads the external JMS known compound database adapter.
ext_Experiment_load_annotation_db <- function(self, src = "hmdb4") {
  factory <- .experiment_get(self, "knownCompoundDatabase_factory")
  if (!is.function(factory)) {
    stop(
      "LC annotation requires a knownCompoundDatabase_factory adapter.",
      call. = FALSE
    )
  }
  self$KCD <- factory(src = src)
  invisible(NULL)
}

# Corresponds to ext_Experiment.db_mass_calibrate: correct all feature m/z according to the system bias estimated by the database.
ext_Experiment_db_mass_calibrate <- function(
    self, required_calibrate_threshold = 0.000002) {
  rows <- vapply(
    self$CMAP$good_reference_landmark_peaks,
    `[[`, numeric(1), "ref_id_num"
  )
  mz_landmarks <- self$CMAP$MassGrid$mz[as.integer(rows) + 1L]
  mass_accuracy_ratio <- .experiment_call(
    self$KCD,
    "evaluate_mass_accuracy_ratio",
    mz_landmarks,
    mode = self$mode,
    mz_tolerance_ppm = 10
  )
  if (!is.null(mass_accuracy_ratio) && length(mass_accuracy_ratio) > 0L &&
      !is.na(mass_accuracy_ratio) && mass_accuracy_ratio != 0) {
    if (abs(mass_accuracy_ratio) > required_calibrate_threshold) {
      message(
        "Mass shift is greater than ",
        sprintf("%2.1f", required_calibrate_threshold * 1e6),
        " ppm. Correction applied."
      )
      correction <- mass_accuracy_ratio + 1
      for (ii in seq_along(self$CMAP$FeatureList)) {
        self$CMAP$FeatureList[[ii]]$mz <-
          self$CMAP$FeatureList[[ii]]$mz / correction
        self$CMAP$FeatureList[[ii]]$mz_corrected_by_division <- correction
      }
    }
  } else {
    message(
      "Mass accuracy check is skipped, too few mz_landmarks (",
      length(mz_landmarks), ") matched."
    )
  }
  invisible(NULL)
}

# Corresponds to ext_Experiment.append_orphans_to_epmCpds: wrap unannotated features into orphan empCpd.
ext_Experiment_append_orphans_to_epmCpds <- function(self, dict_empCpds) {
  all_feature_ids <- unlist(lapply(dict_empCpds, function(compound) {
    vapply(
      compound$MS1_pseudo_Spectra,
      function(peak) as.character(peak$id_number),
      character(1)
    )
  }), use.names = FALSE)
  orphans <- Filter(function(peak) {
    !(as.character(peak$id_number) %in% all_feature_ids)
  }, self$CMAP$FeatureList)
  new_id <- length(dict_empCpds) + 100000L
  for (peak in orphans) {
    dict_empCpds[[as.character(new_id)]] <- list(
      interim_id = new_id,
      neutral_formula_mass = "",
      neutral_formula = "",
      MS1_pseudo_Spectra = list(peak),
      ion_relation = NULL,
      modification = NULL
    )
    new_id <- new_id + 1L
  }
  dict_empCpds
}

# Corresponds to ext_Experiment.export_peak_annotation: Generate feature-by-feature annotation TSV.
ext_Experiment_export_peak_annotation <- function(
    self, dict_empCpds, KCD, export_file_name_prefix) {
  header <- c(
    "[peak]id_number", "mz", "rtime", "apex(scan number)",
    "[EmpCpd]interim_id", "[EmpCpd]ion_relation", "neutral_formula",
    "neutral_formula_mass", "name_1st_guess", "matched_DB_shorts",
    "matched_DB_records"
  )
  lines <- paste(header, collapse = "\t")
  for (compound in dict_empCpds) {
    first_name <- matched_shorts <- matched_records <- ""
    matches <- .experiment_get(compound, "list_matches", list())
    if (length(matches) > 0L) {
      first_match_id <- as.character(matches[[1L]][[1L]])
      db_record <- .experiment_get(KCD, "mass_indexed_compounds")[[first_match_id]]
      first_name <- db_record$compounds[[1L]]$name
      matched_shorts <- paste(vapply(matches, function(match) {
        paste0("(", .experiment_call(
          KCD, "short_report_emp_cpd", match[[1L]]
        ), ")")
      }, character(1)), collapse = ", ")
      matched_records <- paste(vapply(matches, function(match) {
        # Emulates Python str(tuple): character values are enclosed in single quotes, numeric values remain normal text.
        values <- vapply(match, function(value) {
          if (is.character(value)) paste0("'", value, "'") else {
            as.character(value)
          }
        }, character(1))
        paste0("(", paste(values, collapse = ", "), ")")
      }, character(1)), collapse = ", ")
    }
    for (peak in compound$MS1_pseudo_Spectra) {
      row <- c(
        peak$id_number, peak$mz, peak$rtime, peak$apex,
        compound$interim_id, .experiment_get(peak, "ion_relation", ""),
        compound$neutral_formula, compound$neutral_formula_mass,
        first_name, matched_shorts, matched_records
      )
      lines <- c(lines, paste(as.character(row), collapse = "\t"))
    }
  }
  outfile <- file.path(
    self$parameters$outdir, paste0(export_file_name_prefix, ".tsv")
  )
  writeLines(lines, outfile, useBytes = TRUE)
  message(
    "\nAnnotation of ", length(dict_empCpds),
    " Empirical compounds was written to ", outfile, ".\n"
  )
  invisible(NULL)
}

# Corresponds to ext_Experiment.select_unique_compound_features: select a representative peak for each empCpd.
ext_Experiment_select_unique_compound_features <- function(self, dict_empCpds) {
  selected <- list()
  for (compound_id in names(dict_empCpds)) {
    compound <- dict_empCpds[[compound_id]]
    peaks <- compound$MS1_pseudo_Spectra
    if (length(peaks) == 1L) {
      best_peak <- peaks[[1L]]
      relation <- "singleton"
    } else {
      numeric_quality <- all(vapply(peaks, function(peak) {
        is.numeric(peak$peak_area) && is.numeric(peak$goodness_fitting)
      }, logical(1)))
      if (!numeric_quality) {
        best_peak <- peaks[[1L]]
      } else {
        areas <- vapply(peaks, `[[`, numeric(1), "peak_area")
        goodness <- vapply(peaks, `[[`, numeric(1), "goodness_fitting")
        ordering <- order(areas, goodness, decreasing = TRUE)
        best_peak <- peaks[[ordering[[1L]]]]
      }
      relation <- .experiment_get(best_peak, "ion_relation", "")
    }
    selected[[as.character(best_peak$id_number)]] <- list(
      compound$interim_id, compound$neutral_formula, relation
    )
  }
  self$selected_unique_features <- selected
  invisible(NULL)
}

# Corresponds to ext_Experiment.export_feature_tables: export full, preferred and unique tables.
ext_Experiment_export_feature_tables <- function(
    self, `_snr` = 2, `_peak_shape` = 0.7, `_cSelectivity` = 0.7) {
  if (isTRUE(self$parameters$drop_unaligned_samples)) {
    good_samples <- vapply(Filter(function(sample) {
      isTRUE(sample$is_rt_aligned)
    }, self$all_samples), `[[`, character(1), "name")
  } else {
    good_samples <- vapply(self$all_samples, `[[`, character(1), "name")
  }
  self$number_of_samples <- length(good_samples)
  all_sample_names <- vapply(self$all_samples, `[[`, character(1), "name")
  self$dropped_sample_names <- setdiff(all_sample_names, good_samples)

  feature_table <- self$CMAP$FeatureTable
  missing_samples <- setdiff(good_samples, names(feature_table))
  if (length(missing_samples) > 0L) {
    stop("FeatureTable is missing sample columns: ", paste(missing_samples, collapse = ", "))
  }
  detection_counts <- if (length(good_samples) == 0L) {
    integer(nrow(feature_table))
  } else {
    rowSums(as.data.frame(feature_table[good_samples]) > 1, na.rm = TRUE)
  }
  self$CMAP$FeatureTable$detection_counts <- detection_counts
  use_columns <- c(
    "id_number", "mz", "rtime", "rtime_left_base", "rtime_right_base",
    "parent_masstrack_id", "peak_area", "cSelectivity", "goodness_fitting",
    "snr", "detection_counts", good_samples
  )
  filtered <- self$CMAP$FeatureTable[use_columns]
  filtered$mz <- round(filtered$mz, 4)
  filtered$rtime <- round(filtered$rtime, 2)
  filtered$rtime_left_base <- round(filtered$rtime_left_base, 2)
  filtered$rtime_right_base <- round(filtered$rtime_right_base, 2)
  filtered$cSelectivity <- round(filtered$cSelectivity, 2)
  filtered$goodness_fitting <- round(filtered$goodness_fitting, 2)

  full_path <- file.path(
    self$parameters$outdir, "export",
    paste0("full_", self$parameters$output_feature_table)
  )
  utils::write.table(
    filtered, full_path, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = ""
  )
  message(
    "\nFeature table (", nrow(filtered), " x ", length(good_samples),
    ") was written to ", full_path, "."
  )

  preferred <- filtered[
    filtered$detection_counts > 0 & filtered$snr > `_snr` &
      filtered$goodness_fitting > `_peak_shape` &
      filtered$cSelectivity > `_cSelectivity`,
    , drop = FALSE
  ]
  preferred_path <- file.path(
    self$parameters$outdir,
    paste0("preferred_", self$parameters$output_feature_table)
  )
  utils::write.table(
    preferred, preferred_path, sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = ""
  )
  message(
    "\nFiltered Feature table (", nrow(preferred), " x ",
    length(good_samples), ") was written to ", preferred_path, ".\n"
  )

  if (isTRUE(self$parameters$anno) && !identical(self$parameters$workflow, "GC")) {
    selected_ids <- names(self$selected_unique_features)
    unique_table <- preferred[
      as.character(preferred$id_number) %in% selected_ids,
      , drop = FALSE
    ]
    metadata <- lapply(as.character(unique_table$id_number), function(id) {
      self$selected_unique_features[[id]]
    })
    unique_table <- cbind(
      unique_table[seq_len(min(3L, ncol(unique_table)))],
      empCpd = vapply(metadata, `[[`, numeric(1), 1L),
      neutral_formula = vapply(metadata, `[[`, character(1), 2L),
      ion_relation = vapply(metadata, `[[`, character(1), 3L),
      unique_table[if (ncol(unique_table) > 3L) 4:ncol(unique_table) else integer()]
    )
    unique_path <- file.path(
      self$parameters$outdir, "export",
      paste0("unique_compound__", self$parameters$output_feature_table)
    )
    utils::write.table(
      unique_table, unique_path, sep = "\t", row.names = FALSE,
      col.names = TRUE, quote = FALSE, na = ""
    )
    message(
      "Unique compound table (", nrow(unique_table), " x ",
      length(good_samples), ") was written to ", unique_path, ".\n"
    )
  }
  invisible(NULL)
}

# Corresponds to ext_Experiment.export_log: Supplement sample statistics and write project.json.
ext_Experiment_export_log <- function(self) {
  .experiment_set_parameter(self, "number_of_samples", self$number_of_samples)
  .experiment_set_parameter(self, "number_scans", self$number_scans)
  .experiment_set_parameter(
    self, "dropped_samples", .experiment_get(self, "dropped_sample_names", list())
  )
  .experiment_write_json(
    self$parameters,
    file.path(self$parameters$outdir, "project.json")
  )
  invisible(NULL)
}

# Corresponds to ext_Experiment.export_readme: Write Python default description text.
ext_Experiment_export_readme <- function(self) {
  connection <- file(
    file.path(self$parameters$outdir, "README.txt"),
    open = "wb"
  )
  on.exit(close(connection), add = TRUE)
  writeChar(.experiment_readme, connection, eos = NULL, useBytes = TRUE)
  invisible(NULL)
}

# Use a constructor alias close to the Python class name.
ext_Experiment <- ext_Experiment__init__

# Preserves the compatibility constructor used by early smoke_check.R; it is not part of Python's 20 defs.
new_experiment <- function(samples = list(), parameters = default_parameters()) {
  list(samples = samples, parameters = parameters)
}
