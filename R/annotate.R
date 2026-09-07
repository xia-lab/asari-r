# Corresponds to Python asari/annotate.py: LC/GC project annotation entry and result export.

# Call the environment or list adapter method uniformly.
.annotation_call <- function(object, method, ...) {
  callback <- if (is.environment(object)) object[[method]] else object[[method]]
  if (!is.function(callback)) stop("Annotation adapter lacks method: ", method)
  callback(...)
}

# Use jsonlite output and recursively convert the matrix to correspond to NpEncoder.
.annotation_json_ready <- function(value) {
  if (is.matrix(value)) return(lapply(seq_len(nrow(value)), function(ii) as.list(value[ii, ])))
  if (is.environment(value)) return(.annotation_json_ready(as.list(value)))
  if (is.list(value)) return(lapply(value, .annotation_json_ready))
  value
}
.annotation_write_json <- function(value, outfile, pretty = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Writing annotation JSON requires jsonlite.")
  jsonlite::write_json(
    .annotation_json_ready(value), outfile, pretty = pretty,
    auto_unbox = TRUE, na = "null", null = "null"
  )
  invisible(NULL)
}

# Corresponds to annotate_project: Create a timestamp directory and dispatch LC or GC workflow.
annotate_project <- function(infile, parameters) {
  now <- Sys.time()
  stamp <- paste0(vapply(
    c("%m", "%d", "%H", "%M", "%S"),
    function(pattern) as.integer(format(now, pattern)), integer(1)
  ), collapse = "")
  subdir <- paste("annotation", parameters$project_name, stamp, sep = "_")
  outdir <- file.path(parameters$outdir, subdir)
  if (!dir.create(outdir, recursive = FALSE, showWarnings = FALSE)) {
    stop("Cannot create annotation directory: ", outdir)
  }
  cat("Annotation directory: ", outdir, ".\n\n", sep = "")
  if (identical(parameters$workflow, "LC")) {
    annotation <- LCMS_Annotation(parameters)
    annotation$parameters$outdir <- outdir
    annotation$annotate_user_featuretable(infile)
  } else if (identical(parameters$workflow, "GC")) {
    annotate_gcms_full(
      infile, outdir, parameters$kovats, parameters$db,
      project_name_handle = parameters$project_name,
      denovo = parameters$denovo,
      ms2_tolerance_in_ppm = parameters$ms2_tolerance_in_ppm,
      ms2_tolerance_in_da = parameters$ms2_tolerance_in_da,
      ri_tolerance = parameters$ri_tolerance,
      score_cutoff_cosine = parameters$score_cutoff_cosine,
      score_cutoff_entropy = parameters$score_cutoff_entropy,
      corr_cutoff = parameters$corr_cutoff,
      max_ri_delta = parameters$max_ri_delta,
      do_mirror_plot = parameters$do_mirror_plot,
      max_core_features = parameters$max_core_features
    )
  }
  .annotation_write_json(parameters, file.path(outdir, "project.json"))
  invisible(NULL)
}

# Generate safe text that can be used for file names, corresponding to pathvalidate.sanitize_filename.
.annotation_safe_filename <- function(value) {
  value <- gsub("[<>:\"/\\\\|?*]", "_", value)
  sub("[. ]+$", "", value)
}

# Corresponds to annotate_gcms_full: full GC-HRMS target annotation and optional de novo workflow.
annotate_gcms_full <- function(
    infile, outdir, KovatsIndex, database_file,
    project_name_handle = "result", low_peak_filter_factor = 100,
    ms2_tolerance_in_ppm = 5, ms2_tolerance_in_da = 0.005,
    ri_tolerance = 50, cosine_penalty = 1, score_cutoff_cosine = 0.5,
    score_cutoff_entropy = 0.4, corr_cutoff = 0.6, min_ri_delta = 1,
    max_ri_delta = 100, hcl_distance_cut = 1, clustering_step_size = 2000,
    feature_distance_filter = NULL, do_mirror_plot = TRUE,
    max_core_features = 20000, denovo = FALSE) {
  invisible(ms2_tolerance_in_ppm)
  invisible(min_ri_delta)
  invisible(max_ri_delta)
  invisible(feature_distance_filter)
  parsed <- read_features_from_asari_table(paste(readLines(infile, warn = FALSE), collapse = "\n"))
  list_features <- append_kovats_index(
    parsed[[2L]], read_fit_KovatsIndex_rtime(KovatsIndex, sep = "\t", frac = 0.3)
  )
  list_features <- list_features[order(-vapply(list_features, `[[`, 0, "peak_area"))]
  dict_features <- stats::setNames(list_features, vapply(list_features, `[[`, "", "id"))
  list_lib_entries <- reformat_gcms_lib(
    load_gcms_dbfile(database_file), filter_factor = low_peak_filter_factor
  )
  dict_lib_entries <- stats::setNames(list_lib_entries, vapply(list_lib_entries, `[[`, "", "id"))
  cat("Imported ", length(list_lib_entries), " compound library entries.\n", sep = "")

  table <- utils::read.delim(infile, row.names = 1L, check.names = FALSE)
  feature_dataframe <- if (ncol(table) > 10L) table[, 11:ncol(table), drop = FALSE] else table[, 0, drop = FALSE]
  matched <- batch_lib_search_score(
    list_lib_entries, list_features, dict_features, feature_dataframe,
    ms2_tolerance_in_da, ri_tolerance, cosine_penalty, corr_cutoff
  )
  curated <- curate_batch_lib_search_result(
    matched, ms2_tolerance_in_da, score_cutoff_cosine, score_cutoff_entropy
  )
  list_empCpds <- curated[[1L]]
  feature_anno_list <- curated[[2L]]
  write_tsv_feature_anno(
    feature_anno_list, dict_features, dict_lib_entries,
    file.path(outdir, paste0("Features_", project_name_handle, ".tsv"))
  )
  write_tsv_empCpd_anno(
    list_empCpds, dict_features, dict_lib_entries,
    file.path(outdir, paste0("empCpds_", project_name_handle, ".tsv"))
  )
  cat("\nDone targeted annotation.\n")
  unique_features <- unique(vapply(feature_anno_list, `[[`, "", "feature"))
  cat(sprintf(
    "Exported tsv results for %d annotated empCpds and %d unique features.\n\n",
    length(list_empCpds), length(unique_features)
  ))

  if (isTRUE(do_mirror_plot)) {
    cat("Exporting PDF mirror plots..\n\n")
    mirror_dir <- file.path(outdir, "mirrorplots")
    dir.create(mirror_dir, recursive = TRUE, showWarnings = FALSE)
    plotter <- get0("mirror_plot", mode = "function", inherits = TRUE)
    if (!is.function(plotter)) {
      warning("mirror_plot is unavailable; mirror PDF export was skipped.", call. = FALSE)
    } else {
      for (entry in list_empCpds) {
        score <- sprintf(
          "\n(cosine score: %.2f, entropy score: %.2f)",
          entry$cosine_score, entry$entropy_score
        )
        title <- paste0(entry$id, "__", entry$name)
        plotter(
          entry$peaks_as_features,
          dict_lib_entries[[entry$lib_entry_id]]$peaks,
          match_tol = ms2_tolerance_in_da,
          colors = c("blue", "tab:red", "black"),
          title = paste0(title, score),
          outfile = file.path(mirror_dir, paste0(.annotation_safe_filename(title), ".pdf"))
        )
      }
    }
  }

  serialized <- serialize_annotated_empCpds(list_empCpds)
  json_file <- file.path(outdir, paste0(project_name_handle, "_annotated_pseudospectra.json"))
  .annotation_write_json(serialized, json_file)
  json_pseudospectra_to_msp(serialized, sub(".json", ".msp", json_file, fixed = TRUE))
  cat(project_name_handle, "_annotated_pseudospectra was written to JSON and MSP formats.\n\n", sep = "")

  if (isTRUE(denovo)) {
    cat("De novo construction of pseudospectra (deconvolution)...\n")
    core <- unique(vapply(feature_anno_list, `[[`, "", "feature"))
    built <- iterative_build_pseudospectra_by_hcl(
      list_features, feature_dataframe, core, clustering_step_size,
      hcl_distance_cut, ri_tolerance, corr_cutoff, max_core_features
    )
    spectra <- lapply(built[[1L]], port_pseudospectrum_to_json, normalize_intensity = TRUE)
    denovo_file <- file.path(outdir, paste0(project_name_handle, "_denovo_pseudospectra.json"))
    .annotation_write_json(spectra, denovo_file)
    json_pseudospectra_to_msp(spectra, sub(".json", ".msp", denovo_file, fixed = TRUE))
    cat("Exported ", length(spectra), " de novo pseudospectra to JSON and MSP formats.\n\n", sep = "")
  }
  invisible(NULL)
}

# Create a variable environment object corresponding to Python LCMS_Annotation.
LCMS_Annotation <- function(parameters) LCMS_Annotation__init__(parameters)

# Corresponds to LCMS_Annotation.__init__.
LCMS_Annotation__init__ <- function(parameters) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("LCMS_Annotation", "environment")
  self$parameters <- parameters
  self$list_features <- list()
  self$KCD <- NULL
  self$annotate_user_featuretable <- function(infile) LCMS_Annotation_annotate_user_featuretable(self, infile)
  self$load_annotation_db <- function(src = "hmdb4") LCMS_Annotation_load_annotation_db(self, src)
  self$db_mass_calibrate <- function(max_features = 1000, required_calibrate_threshold = 0.000002) {
    LCMS_Annotation_db_mass_calibrate(self, max_features, required_calibrate_threshold)
  }
  self$append_orphans_to_epmCpds <- function(dict_empCpds) LCMS_Annotation_append_orphans_to_epmCpds(self, dict_empCpds)
  self$export_peak_annotation <- function(dict_empCpds, KCD, prefix) {
    LCMS_Annotation_export_peak_annotation(self, dict_empCpds, KCD, prefix)
  }
  self$select_unique_compound_features <- function(dict_empCpds) {
    LCMS_Annotation_select_unique_compound_features(self, dict_empCpds)
  }
  self
}

# Corresponds to the LCMS_Annotation.annotate_user_featuretable three-stage JMS annotation process.
LCMS_Annotation_annotate_user_featuretable <- function(self, infile) {
  parsed <- read_features_from_asari_table(paste(readLines(infile, warn = FALSE), collapse = "\n"))
  features <- lapply(parsed[[2L]], function(feature) {
    feature$representative_intensity <- feature$peak_area
    feature
  })
  self$list_features <- features[order(-vapply(features, `[[`, 0, "peak_area"))]
  self$load_annotation_db()
  self$db_mass_calibrate()
  factory <- getOption("asariR.ExperimentalEcpdDatabase_factory")
  if (!is.function(factory)) stop("LC annotation requires ExperimentalEcpdDatabase_factory adapter.")
  EED <- factory(
    mode = self$parameters$mode,
    mz_tolerance_ppm = self$parameters$mz_tolerance_ppm,
    rt_tolerance = self$parameters$khipu_rtime_tolerance
  )
  EED$adduct_patterns <- if (identical(self$parameters$mode, "pos")) {
    adduct_search_patterns_pos
  } else adduct_search_patterns_neg
  EED$isotope_search_patterns <- isotope_search_patterns
  EED$extended_adducts <- extended_adducts
  .annotation_call(EED, "build_from_list_peaks", self$list_features)
  .annotation_call(EED, "extend_empCpd_annotation", self$KCD)
  .annotation_call(EED, "annotate_singletons", self$KCD)
  EED$dict_empCpds <- self$append_orphans_to_epmCpds(EED$dict_empCpds)
  self$export_peak_annotation(EED$dict_empCpds, self$KCD, "Feature_annotation")
  outfile <- file.path(self$parameters$outdir, "Annotated_empiricalCompounds.json")
  .annotation_write_json(EED$dict_empCpds, outfile)
  cat("JSON version was written to ", outfile, ".\n\n", sep = "")
  invisible(NULL)
}

# Corresponds to LCMS_Annotation.load_annotation_db: JMS database provided by an injectable factory.
LCMS_Annotation_load_annotation_db <- function(self, src = "hmdb4") {
  factory <- getOption("asariR.knownCompoundDatabase_factory")
  if (!is.function(factory)) {
    stop(
      "LC annotation requires knownCompoundDatabase_factory; Python pickle data alone do not implement JMS methods.",
      call. = FALSE
    )
  }
  self$KCD <- factory(src = src)
  invisible(NULL)
}

# The first db_mass_calibrate definition in the corresponding class only has a docstring; this definition is subsequently overridden by Python.
LCMS_Annotation_db_mass_calibrate_documented <- function(
    self, max_features = 1000, required_calibrate_threshold = 0.000002) {
  invisible(self)
  invisible(max_features)
  invisible(required_calibrate_threshold)
  invisible(NULL)
}

# Corresponds to the final effective LCMS_Annotation.db_mass_calibrate.
LCMS_Annotation_db_mass_calibrate <- function(
    self, max_features = 1000, required_calibrate_threshold = 0.000002) {
  landmarks <- utils::head(vapply(self$list_features, `[[`, 0, "mz"), max_features)
  ratio <- .annotation_call(
    self$KCD, "evaluate_mass_accuracy_ratio", landmarks,
    mode = self$parameters$mode, mz_tolerance_ppm = 10
  )
  if (!is.null(ratio) && length(ratio) && !is.na(ratio) && ratio != 0) {
    if (abs(ratio) > required_calibrate_threshold) {
      cat(sprintf(
        "Mass shift is greater than %2.1f ppm. Correction applied.\n",
        required_calibrate_threshold * 1e6
      ))
      correction <- ratio + 1
      self$list_features <- lapply(self$list_features, function(feature) {
        feature$mz <- feature$mz / correction
        feature$mz_corrected_by_division <- correction
        feature
      })
    }
  } else {
    cat("Mass accuracy check is skipped, too few mz_landmarks (", length(landmarks), ") matched.\n", sep = "")
  }
  invisible(NULL)
}

# Corresponds to append_orphans_to_epmCpds.
LCMS_Annotation_append_orphans_to_epmCpds <- function(self, dict_empCpds) {
  assigned <- unlist(lapply(dict_empCpds, function(compound) {
    vapply(compound$MS1_pseudo_Spectra, function(peak) as.character(peak$id_number), "")
  }), use.names = FALSE)
  orphans <- Filter(function(peak) !(as.character(peak$id_number) %in% assigned), self$list_features)
  new_id <- length(dict_empCpds) + 100000L
  for (peak in orphans) {
    dict_empCpds[[as.character(new_id)]] <- list(
      interim_id = new_id, neutral_formula_mass = "", neutral_formula = "",
      MS1_pseudo_Spectra = list(peak), ion_relation = NULL, modification = NULL
    )
    new_id <- new_id + 1L
  }
  dict_empCpds
}

# Corresponds to export_peak_annotation: output the peak-by-peak JMS annotation table.
LCMS_Annotation_export_peak_annotation <- function(self, dict_empCpds, KCD, export_file_name_prefix) {
  header <- c(
    "[peak]id_number", "mz", "rtime", "apex(scan number)",
    "[EmpCpd]interim_id", "[EmpCpd]ion_relation", "neutral_formula",
    "neutral_formula_mass", "name_1st_guess", "matched_DB_shorts", "matched_DB_records"
  )
  lines <- paste(header, collapse = "\t")
  for (compound in dict_empCpds) {
    first_name <- shorts <- records <- ""
    matches <- compound$list_matches %||% list()
    if (length(matches)) {
      first_id <- as.character(matches[[1L]][[1L]])
      first_name <- KCD$mass_indexed_compounds[[first_id]]$compounds[[1L]]$name
      shorts <- paste(vapply(matches, function(match) {
        paste0("(", .annotation_call(KCD, "short_report_emp_cpd", match[[1L]]), ")")
      }, ""), collapse = ", ")
      records <- paste(vapply(matches, function(match) paste0("(", paste(match, collapse = ", "), ")"), ""), collapse = ", ")
    }
    for (peak in compound$MS1_pseudo_Spectra) {
      row <- c(
        peak$id_number, peak$mz, peak$rtime, peak$apex, compound$interim_id,
        peak$ion_relation %||% "", compound$neutral_formula,
        compound$neutral_formula_mass, first_name, shorts, records
      )
      lines <- c(lines, paste(row, collapse = "\t"))
    }
  }
  outfile <- file.path(self$parameters$outdir, paste0(export_file_name_prefix, ".tsv"))
  writeLines(lines, outfile, useBytes = TRUE)
  cat(sprintf(
    "\nAnnotation of %d Empirical compounds was written to %s.\n\n",
    length(dict_empCpds), outfile
  ))
  invisible(NULL)
}

# Corresponds to select_unique_compound_features: select the highest peak of comprehensive area for each empCpd.
LCMS_Annotation_select_unique_compound_features <- function(self, dict_empCpds) {
  selected <- list()
  for (interim_id in names(dict_empCpds)) {
    compound <- dict_empCpds[[interim_id]]
    peaks <- compound$MS1_pseudo_Spectra
    if (length(peaks) == 1L) {
      selected[[as.character(peaks[[1L]]$id_number)]] <- list(
        interim_id, compound$neutral_formula, "singleton"
      )
    } else {
      area <- vapply(peaks, `[[`, 0, "peak_area")
      shape <- vapply(peaks, `[[`, 0, "goodness_fitting")
      ties <- which(area == max(area) & shape == max(shape[area == max(area)]))
      best <- peaks[[ties[[1L]]]]
      selected[[as.character(best$id_number)]] <- list(
        interim_id, compound$neutral_formula, best$ion_relation %||% ""
      )
    }
  }
  self$selected_unique_features <- selected
  invisible(NULL)
}

# Corresponding module level annotate_user_featuretable: independent khipu pre-annotation entry.
annotate_user_featuretable <- function(infile, parameters) {
  mode <- parameters$mode
  patterns <- if (identical(mode, "pos")) adduct_search_patterns_pos else adduct_search_patterns_neg
  features <- read_features_from_asari_table(
    paste(readLines(infile, warn = FALSE), collapse = "\n")
  )[[2L]]
  features <- lapply(features, function(feature) {
    feature$representative_intensity <- feature$peak_area
    feature
  })
  converter <- getOption("asariR.peaklist_to_khipu_list")
  if (!is.function(converter)) stop("Standalone khipu annotation requires peaklist_to_khipu_list adapter.")
  converted <- converter(
    features, isotope_search_patterns = isotope_search_patterns,
    adduct_search_patterns = patterns, extended_adducts = extended_adducts,
    mz_tolerance_ppm = parameters$mz_tolerance_ppm,
    rt_tolerance = parameters$khipu_rtime_tolerance,
    mode = mode, charges = c(1L, 2L, 3L)
  )
  outfile <- file.path(
    parameters$outdir,
    paste0(parameters$project_name, "_PreAnnotated_list_khipus.json")
  )
  .annotation_write_json(converted[[1L]], outfile)
  invisible(NULL)
}

# Corresponds to annoate_by_standards with the spelling unchanged.
annoate_by_standards <- function(list_features, cpdLib) {
  matched <- list_match_lcms_features(list_features, cpdLib, mz_ppm = 5, rt_tolerance = 30)
  dict_f <- stats::setNames(list_features, vapply(list_features, `[[`, "", "id"))
  dict_lib <- stats::setNames(cpdLib, vapply(cpdLib, `[[`, "", "id"))
  invisible(dict_f)
  invisible(dict_lib)
  matched
}

# Corresponds to get_concise_annotation: each class library only retains the first match.
get_concise_annotation <- function(feature_id, matched_list, feature_dict, lib_dict) {
  feature <- feature_dict[[feature_id]]
  result <- list(id = feature_id, mz = feature$mz, rtime = feature$rtime)
  authentic <- matched_list[startsWith(matched_list, "v2r2024")]
  csm <- matched_list[startsWith(matched_list, "r1_")]
  if (length(authentic)) {
    id <- authentic[[1L]]
    info <- lib_dict[[id]]
    result <- c(result, list(
      lib_id = id, lib_name = info$name, lib_mz = info$mz,
      lib_rtime = info$rtime, lib_identifier = info$identifier,
      lib_ion = info$ion, lib_isotope = info$isotope
    ))
  } else {
    result <- c(result, stats::setNames(as.list(rep("", 7L)), c(
      "lib_id", "lib_name", "lib_mz", "lib_rtime", "lib_identifier", "lib_ion", "lib_isotope"
    )))
  }
  if (length(csm)) {
    id <- csm[[1L]]
    info <- lib_dict[[id]]
    result <- c(result, list(
      CSMF_ID = id, CSM_ion = info$ion_csm,
      CSM_top_recommendation_name = info$top_recommendation_name,
      CSM_top_recommendation_score = info$top_recommendation_score,
      CSM_HMDB = info$HMDB
    ))
  } else {
    result <- c(result, stats::setNames(as.list(rep("", 5L)), c(
      "CSMF_ID", "CSM_ion", "CSM_top_recommendation_name",
      "CSM_top_recommendation_score", "CSM_HMDB"
    )))
  }
  result
}

# Corresponds to export_combined_anno_table.
export_combined_anno_table <- function(matched, feature_dict, lib_dict, outfile) {
  header <- c(
    "id", "mz", "rtime", "lib_id", "lib_name", "lib_mz", "lib_rtime",
    "lib_identifier", "lib_ion", "lib_isotope", "CSMF_ID", "CSM_ion",
    "CSM_top_recommendation_name", "CSM_top_recommendation_score", "CSM_HMDB"
  )
  lines <- paste(header, collapse = "\t")
  for (feature_id in names(matched)) {
    concise <- get_concise_annotation(feature_id, matched[[feature_id]], feature_dict, lib_dict)
    lines <- c(lines, paste(vapply(header, function(field) as.character(concise[[field]]), ""), collapse = "\t"))
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to export_combined_anno_json.
export_combined_anno_json <- function(matched, feature_dict, lib_dict, outfile) {
  output <- lapply(names(matched), function(feature_id) {
    stats::setNames(list(list(
      Feature = feature_dict[[feature_id]],
      Matched_libs = lapply(matched[[feature_id]], function(id) lib_dict[[id]])
    )), feature_id)
  })
  .annotation_write_json(output, outfile)
  invisible(NULL)
}
