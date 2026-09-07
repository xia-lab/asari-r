# Public annotation and MS2 entry points: all inputs must be explicit objects or file paths.

# Find the actual column name used by the database from a set of possible field names.
.asari_annotation_column <- function(table, candidates) {
  actual <- names(table)
  matched <- match(tolower(candidates), tolower(actual), nomatch = 0L)
  matched <- matched[matched > 0L]
  if (length(matched)) actual[[matched[[1L]]]] else NULL
}

# Safely retrieve the first known field in a list record.
.asari_annotation_value <- function(record, candidates, fallback = NULL) {
  actual <- names(record)
  if (is.null(actual)) return(fallback)
  matched <- match(tolower(candidates), tolower(actual), nomatch = 0L)
  matched <- matched[matched > 0L]
  if (length(matched)) record[[actual[[matched[[1L]]]]]] else fallback
}

# Unify the tabular LC database into the 8 columns required to expose annotation functions.
.asari_normalize_lc_table <- function(table) {
  if (!is.data.frame(table)) stop("LC database must be a data.frame or supported file.", call. = FALSE)
  n <- nrow(table)
  get_text <- function(candidates, fallback = rep("", n)) {
    column <- .asari_annotation_column(table, candidates)
    if (is.null(column)) fallback else as.character(table[[column]])
  }
  get_number <- function(candidates) {
    column <- .asari_annotation_column(table, candidates)
    if (is.null(column)) rep(NA_real_, n) else suppressWarnings(as.numeric(table[[column]]))
  }
  result <- data.frame(
    compound_id = get_text(c("compound_id", "primary_id", "id", "accession")),
    name = get_text(c("name", "compound_name", "title")),
    mz = get_number(c("mz", "precursor_mz", "ion_mz")),
    rtime = get_number(c("rtime", "retention_time", "rt")),
    neutral_mass = get_number(c(
      "neutral_mass", "neutral_formula_mass", "exact_mass", "exactmass", "mw"
    )),
    formula = get_text(c("formula", "neutral_formula", "molecular_formula")),
    inchikey = get_text(c("inchikey", "inchi_key")),
    source = get_text(c("source", "primary_db", "database")),
    stringsAsFactors = FALSE
  )
  # Use reproducible row numbers instead of random values when the database does not have stable IDs.
  empty_id <- is.na(result$compound_id) | !nzchar(result$compound_id)
  result$compound_id[empty_id] <- paste0("DB", which(empty_id))
  if (all(is.na(result$mz)) && all(is.na(result$neutral_mass))) {
    stop("LC database requires an mz or neutral_mass/exact_mass column.", call. = FALSE)
  }
  result
}

# Flatten asari's mass_indexed_compounds.pickle into a normal data.frame.
.asari_flatten_lc_pickle <- function(database) {
  rows <- list()
  counter <- 0L
  for (entry in unname(database)) {
    compounds <- .asari_annotation_value(entry, "compounds", list())
    # Individual mass records may not have compounds, but still retain masses and formulas.
    if (length(compounds) == 0L) compounds <- list(entry)
    for (compound in compounds) {
      counter <- counter + 1L
      rows[[counter]] <- data.frame(
        compound_id = as.character(.asari_annotation_value(
          compound, c("primary_id", "id"), paste0("DB", counter)
        )),
        name = as.character(.asari_annotation_value(compound, "name", "")),
        mz = NA_real_,
        rtime = NA_real_,
        neutral_mass = suppressWarnings(as.numeric(.asari_annotation_value(
          compound, c("neutral_formula_mass", "exact_mass"),
          .asari_annotation_value(entry, "neutral_formula_mass", NA_real_)
        ))),
        formula = as.character(.asari_annotation_value(
          compound, "neutral_formula",
          .asari_annotation_value(entry, "neutral_formula", "")
        )),
        inchikey = as.character(.asari_annotation_value(compound, "inchikey", "")),
        source = as.character(.asari_annotation_value(compound, "primary_db", "HMDB")),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(data.frame(
    compound_id = character(), name = character(), mz = numeric(), rtime = numeric(),
    neutral_mass = numeric(), formula = character(), inchikey = character(),
    source = character(), stringsAsFactors = FALSE
  ))
  do.call(rbind, rows)
}

#' Read LC-MS annotation database
#'
#' Supports data.frame, TSV, CSV, JSON and asari raw Python pickle. Use
#' When `database = "hmdb4"` is used, an available Python decoding pickle is required for the first time.
#' The RDS cache will then be read.
#'
#' @param database Database data.frame, explicit file path or `"hmdb4"`.
#' @param cache Whether to cache an R-formatted copy of the built-in HMDB.
#' @return data.frame with unified fields.
#' @export
asari_load_lc_database <- function(database = "hmdb4", cache = TRUE) {
  if (is.data.frame(database)) return(.asari_normalize_lc_table(database))
  if (!is.character(database) || length(database) != 1L || is.na(database)) {
    stop("database must be a data.frame, one file path, or 'hmdb4'.", call. = FALSE)
  }
  if (tolower(database) %in% c("hmdb4", "hmdb")) {
    cache_dir <- tools::R_user_dir("asariR", "cache")
    cache_file <- file.path(cache_dir, "hmdb4_lc_database.rds")
    if (isTRUE(cache) && file.exists(cache_file)) return(readRDS(cache_file))
    # Find Python asari data resources only if the user explicitly selects hmdb4.
    raw <- load_asari_db_resource("mass_indexed_compounds.pickle")
    result <- .asari_flatten_lc_pickle(raw)
    if (isTRUE(cache)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, cache_file)
    }
    return(result)
  }
  path <- path.expand(database)
  if (!file.exists(path) || dir.exists(path)) stop("Database file does not exist: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  raw <- if (extension %in% c("tsv", "txt")) {
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (extension == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (extension == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading JSON requires jsonlite.")
    jsonlite::fromJSON(path, simplifyDataFrame = TRUE)
  } else if (extension %in% c("pickle", "pkl")) {
    .samples_read_pickle(path)
  } else {
    stop("Supported LC database formats are TSV, CSV, JSON and pickle.", call. = FALSE)
  }
  if (is.data.frame(raw)) .asari_normalize_lc_table(raw) else .asari_flatten_lc_pickle(raw)
}

# Generate a list of theoretical ions for neutral masses consistent with asari's default rules.
.asari_lc_theoretical_ions <- function(database, mode, include_adducts) {
  records <- list()
  counter <- 0L
  for (ii in seq_len(nrow(database))) {
    if (!is.na(database$mz[[ii]])) {
      counter <- counter + 1L
      records[[counter]] <- data.frame(
        db_row = ii, theoretical_mz = database$mz[[ii]], ion_relation = "database_mz"
      )
      next
    }
    neutral <- database$neutral_mass[[ii]]
    if (is.na(neutral)) next
    base <- neutral + if (mode == "pos") PROTON else -PROTON
    labels <- if (mode == "pos") "[M+H]+" else "[M-H]-"
    masses <- base
    if (isTRUE(include_adducts)) {
      patterns <- if (mode == "pos") adduct_search_patterns_pos else adduct_search_patterns_neg
      masses <- c(masses, base + vapply(patterns, function(x) as.numeric(x[[1L]]), 0))
      labels <- c(labels, vapply(patterns, function(x) as.character(x[[2L]]), ""))
    }
    for (jj in seq_along(masses)) {
      counter <- counter + 1L
      records[[counter]] <- data.frame(
        db_row = ii, theoretical_mz = masses[[jj]], ion_relation = labels[[jj]]
      )
    }
  }
  if (length(records) == 0L) {
    return(data.frame(db_row = integer(), theoretical_mz = numeric(), ion_relation = character()))
  }
  do.call(rbind, records)
}

#' Annotate LC-MS features using mass or standards databases
#'
#' This is a candidate annotation native to R that can be run independently: it matches m/z when RT is provided by the database
#' Then match RT and can extend asari default adduct. It does not claim to replace JMS/Khipu
#' Empirical compound groupings.
#'
#' @param feature_table asari feature table data.frame or explicit TSV path.
#' @param database The database passed to [asari_load_lc_database()].
#' @param output Optional explicit TSV output path.
#' @param mode `"pos"` or `"neg"`.
#' @param ppm m/z tolerance.
#' @param rt_tolerance Optional RT tolerance, units consistent with input table.
#' @param include_adducts Whether to extend the default adduct for neutral mass.
#' @return One feature-candidate compound match per row.
#' @export
asari_annotate_lc <- function(
    feature_table,
    database = "hmdb4",
    output = NULL,
    mode = c("pos", "neg"),
    ppm = 5,
    rt_tolerance = NULL,
    include_adducts = TRUE) {
  mode <- match.arg(mode)
  ppm <- .asari_public_positive_number(ppm, "ppm")
  if (!is.null(rt_tolerance)) {
    rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  }
  features <- .asari_public_feature_table(feature_table)
  .asari_public_require_columns(features, c("mz", "rtime"))
  if (!("id_number" %in% names(features))) features$id_number <- paste0("F", seq_len(nrow(features)) - 1L)
  db <- asari_load_lc_database(database)
  ions <- .asari_lc_theoretical_ions(db, mode, include_adducts)
  ions <- ions[order(ions$theoretical_mz), , drop = FALSE]

  # Search ppm windows using ordered mass indexes to avoid full Cartesian products for large HMDBs.
  matches <- list()
  counter <- 0L
  ion_mz <- ions$theoretical_mz
  for (ii in seq_len(nrow(features))) {
    query <- suppressWarnings(as.numeric(features$mz[[ii]]))
    if (is.na(query) || query <= 0 || length(ion_mz) == 0L) next
    delta <- query * ppm * 1e-6
    left <- findInterval(query - delta, ion_mz) + 1L
    right <- findInterval(query + delta, ion_mz)
    if (left > right || left > length(ion_mz) || right < 1L) next
    candidates <- seq.int(max(1L, left), min(length(ion_mz), right))
    for (jj in candidates) {
      db_row <- ions$db_row[[jj]]
      db_rt <- db$rtime[[db_row]]
      rt_delta <- if (is.na(db_rt)) NA_real_ else as.numeric(features$rtime[[ii]]) - db_rt
      if (!is.null(rt_tolerance) && !is.na(rt_delta) && abs(rt_delta) > rt_tolerance) next
      counter <- counter + 1L
      matches[[counter]] <- data.frame(
        feature_id = as.character(features$id_number[[ii]]),
        feature_mz = query,
        feature_rtime = as.numeric(features$rtime[[ii]]),
        compound_id = db$compound_id[[db_row]],
        compound_name = db$name[[db_row]],
        formula = db$formula[[db_row]],
        inchikey = db$inchikey[[db_row]],
        source = db$source[[db_row]],
        ion_relation = ions$ion_relation[[jj]],
        theoretical_mz = ions$theoretical_mz[[jj]],
        ppm_error = (query - ions$theoretical_mz[[jj]]) / ions$theoretical_mz[[jj]] * 1e6,
        database_rtime = db_rt,
        rt_delta = rt_delta,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- if (length(matches)) do.call(rbind, matches) else data.frame(
    feature_id = character(), feature_mz = numeric(), feature_rtime = numeric(),
    compound_id = character(), compound_name = character(), formula = character(),
    inchikey = character(), source = character(), ion_relation = character(),
    theoretical_mz = numeric(), ppm_error = numeric(), database_rtime = numeric(),
    rt_delta = numeric(), stringsAsFactors = FALSE
  )
  if (nrow(result)) result <- result[order(result$feature_id, abs(result$ppm_error)), , drop = FALSE]
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    attr(result, "output_file") <- outfile
  }
  class(result) <- c("asari_annotation", class(result))
  result
}

#' Annotate GC-MS feature table using Kovats index and EI spectral library
#'
#' @param feature_table an explicit asari full feature table TSV.
#' @param kovats Explicit Kovats calibrated TSV.
#' @param library Explicit MSP or JSON EI spectral library.
#' @param output_dir Comment output directory.
#' @param project_name Output filename identifier.
#' @param denovo Whether to construct de novo pseudo-spectrum at the same time.
#' @param mirror_plots Whether to output mirrored plots; requires an available drawing backend.
#' @param... High-level parameters passed to the underlying GC annotation.
#' @return A named list containing annotation output paths.
#' @export
asari_annotate_gc <- function(
    feature_table,
    kovats,
    library,
    output_dir,
    project_name = "result",
    denovo = FALSE,
    mirror_plots = FALSE,
    ...) {
  required <- c(feature_table = feature_table, kovats = kovats, library = library)
  normalized <- vapply(required, function(path) {
    path <- path.expand(path)
    if (!file.exists(path) || dir.exists(path)) stop("Required GC file does not exist: ", path, call. = FALSE)
    normalizePath(path, mustWork = TRUE)
  }, "")
  output_dir <- path.expand(output_dir)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create GC annotation directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  # The optional mirror image is turned off by default in the public portal, ensuring that the export of tables and spectra is not blocked by the drawing package.
  annotate_gcms_full(
    infile = normalized[["feature_table"]], outdir = output_dir,
    KovatsIndex = normalized[["kovats"]], database_file = normalized[["library"]],
    project_name_handle = project_name, denovo = denovo,
    do_mirror_plot = mirror_plots, ...
  )
  result <- list(
    output_dir = output_dir,
    feature_annotations = file.path(output_dir, paste0("Features_", project_name, ".tsv")),
    empirical_compounds = file.path(output_dir, paste0("empCpds_", project_name, ".tsv")),
    annotated_json = file.path(output_dir, paste0(project_name, "_annotated_pseudospectra.json")),
    annotated_msp = file.path(output_dir, paste0(project_name, "_annotated_pseudospectra.msp"))
  )
  if (isTRUE(denovo)) {
    result$denovo_json <- file.path(output_dir, paste0(project_name, "_denovo_pseudospectra.json"))
    result$denovo_msp <- file.path(output_dir, paste0(project_name, "_denovo_pseudospectra.msp"))
  }
  class(result) <- c("asari_gc_annotation", "list")
  result
}

#' Extract MS1, MS2 and higher level spectra from mzML
#'
#' @param input An explicit mzML file.
#' @param output Optional JSON output path.
#' @param min_intensity Minimum peak intensity.
#' @param ms2_peak_limit The maximum number of peaks retained in each MS2 spectrum; `NULL` means no limit.
#' @return List containing `ms1`, `ms2` and `other`.
#' @export
asari_extract_ms2 <- function(input, output = NULL, min_intensity = 1000, ms2_peak_limit = 50L) {
  infile <- .asari_public_one_mzml(input)
  min_intensity <- .asari_public_positive_number(min_intensity, "min_intensity")
  if (!is.null(ms2_peak_limit)) {
    ms2_peak_limit <- as.integer(.asari_public_positive_number(ms2_peak_limit, "ms2_peak_limit"))
  }
  extracted <- extract_all_spectra_form_file(infile, min_intensity, ms2_peak_limit)
  names(extracted) <- c("ms1", "ms2", "other")
  extracted$input_file <- infile
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    jsonlite::write_json(extracted[c("ms1", "ms2", "other")], outfile, pretty = TRUE, auto_unbox = TRUE)
    extracted$output_file <- outfile
  }
  class(extracted) <- c("asari_ms2_extract", "list")
  extracted
}

#' Match MS2 spectra to MS1 features by m/z and RT
#'
#' @param feature_table Explicit asari full feature table TSV.
#' @param ms2_files One or more explicit mzML files.
#' @param output Explicit JSON output path.
#' @param rt_tolerance RT tolerance, unit is seconds.
#' @param ppm m/z tolerance.
#' @return List of best MS2 spectra organized by MS1 feature ID.
#' @export
asari_match_ms2 <- function(feature_table, ms2_files, output, rt_tolerance = 30, ppm = 5) {
  feature_table <- path.expand(feature_table)
  if (!file.exists(feature_table) || dir.exists(feature_table)) stop("Feature table does not exist: ", feature_table)
  if (!is.character(ms2_files) || length(ms2_files) == 0L) stop("ms2_files must contain explicit mzML paths.")
  files <- vapply(ms2_files, function(path) .asari_public_one_mzml(path), "")
  outfile <- .asari_public_output_file(output)
  match_ms2files_to_features(
    normalizePath(feature_table, mustWork = TRUE), files,
    rt_tol = .asari_public_positive_number(rt_tolerance, "rt_tolerance"),
    ppm_tol = .asari_public_positive_number(ppm, "ppm"), output_json = outfile
  )
}

# Unify MSP, MGF or JSON spectral library into ID, name, precursor m/z and two peak matrices.
.asari_spectral_library <- function(library) {
  if (is.list(library) && !is.character(library)) raw <- library else {
    path <- path.expand(library)
    if (!file.exists(path) || dir.exists(path)) stop("Spectral library does not exist: ", path)
    extension <- tolower(tools::file_ext(path))
    raw <- if (extension == "msp") parse_msp(path) else if (extension == "mgf") {
      parse_mgf(path)
    } else if (extension == "json") {
      jsonlite::fromJSON(path, simplifyVector = FALSE)
    } else stop("Spectral library must be MSP, MGF or JSON.")
  }
  lapply(seq_along(raw), function(ii) {
    entry <- raw[[ii]]
    params <- .asari_annotation_value(entry, "params", list())
    pepmass <- .asari_annotation_value(params, "PEPMASS", NULL)
    precursor <- .asari_annotation_value(
      entry, c("precursor_mz", "precursormz", "precursor m/z"),
      if (!is.null(pepmass)) pepmass[[1L]] else NA_real_
    )
    if (is.na(suppressWarnings(as.numeric(precursor)))) {
      precursor <- .asari_annotation_value(entry, c("exactmass", "exact_mass", "mw"), NA_real_)
    }
    peaks <- if (!is.null(entry$mzs) && !is.null(entry$intensities)) {
      cbind(as.numeric(entry$mzs), as.numeric(entry$intensities))
    } else if (length(entry$peaks) && is.list(entry$peaks[[1L]]) &&
               !is.null(entry$peaks[[1L]]$mz)) {
      do.call(rbind, lapply(entry$peaks, function(peak) {
        c(as.numeric(peak$mz), as.numeric(peak$intensity))
      }))
    } else .ms2_peak_matrix(entry$peaks)
    list(
      id = as.character(.asari_annotation_value(entry, c("id", "accession"), paste0("LIB", ii))),
      name = as.character(.asari_annotation_value(
        entry, c("name", "title"), .asari_annotation_value(params, "TITLE", paste0("LIB", ii))
      )),
      precursor_mz = suppressWarnings(as.numeric(precursor)),
      peaks = peaks,
      metadata = entry
    )
  })
}

#' Search MS2 spectra using local spectral library
#'
#' @param spectra mzML file, [asari_extract_ms2()] result or MS2 spectrum list.
#' @param library MSP, MGF, JSON path or parsed spectral library list.
#' @param output Optional TSV output path.
#' @param ms1_tolerance The absolute tolerance of precursor ions is large.
#' @param ms2_tolerance The absolute tolerance of fragment peaks is large.
#' @param min_similarity Minimum similarity.
#' @param method `"entropy"` or `"cosine"`. When there is no external backend, entropy uses local JS similarity implementation.
#' @param min_intensity The lowest peak intensity when extracted from mzML.
#' @return The best spectral library matching data.frame for each query spectrum.
#' @export
asari_search_ms2 <- function(
    spectra,
    library,
    output = NULL,
    ms1_tolerance = 0.01,
    ms2_tolerance = 0.02,
    min_similarity = 0.6,
    method = c("entropy", "cosine"),
    min_intensity = 1000) {
  method <- match.arg(method)
  ms1_tolerance <- .asari_public_positive_number(ms1_tolerance, "ms1_tolerance")
  ms2_tolerance <- .asari_public_positive_number(ms2_tolerance, "ms2_tolerance")
  if (!is.numeric(min_similarity) || length(min_similarity) != 1L || is.na(min_similarity) ||
      min_similarity < 0 || min_similarity > 1) stop("min_similarity must be between 0 and 1.")
  queries <- if (is.character(spectra)) {
    extract_ms2_form_file(.asari_public_one_mzml(spectra), min_intensity)
  } else if (inherits(spectra, "asari_ms2_extract")) {
    spectra$ms2
  } else if (is.list(spectra)) spectra else stop("spectra must be mzML or an MS2 spectrum list.")
  lib <- .asari_spectral_library(library)
  rows <- list()
  for (ii in seq_along(queries)) {
    query <- queries[[ii]]
    precursor <- as.numeric(query$precursor_mz)
    if (length(precursor) == 0L || is.na(precursor)) next
    candidates <- which(vapply(lib, function(entry) {
      !is.na(entry$precursor_mz) && abs(entry$precursor_mz - precursor) <= ms1_tolerance
    }, FALSE))
    if (length(candidates) == 0L) next
    query_peaks <- .ms2_peak_matrix(query$peaks)
    scores <- lapply(candidates, function(jj) {
      reference <- lib[[jj]]$peaks
      if (method == "entropy") {
        list(.gcms_entropy_similarity(query_peaks, reference, ms2_tolerance),
             length(find_matches(reference[, 1L], query_peaks[, 1L], ms2_tolerance)[[1L]]))
      } else cosine_similarity(query_peaks, reference, tolerance = ms2_tolerance)
    })
    similarity <- vapply(scores, function(score) as.numeric(score[[1L]]), 0)
    best_position <- which.max(similarity)
    if (similarity[[best_position]] < min_similarity) next
    best <- candidates[[best_position]]
    rows[[length(rows) + 1L]] <- data.frame(
      query_id = as.character(query$id %||% paste0("sp", ii)),
      query_precursor_mz = precursor,
      query_rtime = as.numeric(query$rtime %||% NA_real_),
      library_id = lib[[best]]$id,
      library_name = lib[[best]]$name,
      library_precursor_mz = lib[[best]]$precursor_mz,
      precursor_delta = precursor - lib[[best]]$precursor_mz,
      method = if (method == "entropy") "local_js_entropy" else "cosine",
      similarity = similarity[[best_position]],
      matched_peaks = as.integer(scores[[best_position]][[2L]]),
      stringsAsFactors = FALSE
    )
  }
  result <- if (length(rows)) do.call(rbind, rows) else data.frame(
    query_id = character(), query_precursor_mz = numeric(), query_rtime = numeric(),
    library_id = character(), library_name = character(), library_precursor_mz = numeric(),
    precursor_delta = numeric(), method = character(), similarity = numeric(),
    matched_peaks = integer(), stringsAsFactors = FALSE
  )
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    attr(result, "output_file") <- outfile
  }
  result
}
