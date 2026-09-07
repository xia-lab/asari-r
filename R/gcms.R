# Corresponds to Python asari/gcms.py: GC-HRMS pseudo-spectrum construction, spectral library retrieval and annotation export.

# R counterpart of NamedTuple PseudoSpectrum.
PseudoSpectrum <- function(id, rtime, RI, rounded_mzs, num_features, members, peaks, annotation) {
  structure(list(
    id = id, rtime = rtime, RI = RI, rounded_mzs = rounded_mzs,
    num_features = as.integer(num_features), members = members,
    peaks = as.matrix(peaks), annotation = annotation
  ), class = c("PseudoSpectrum", "list"))
}

# NamedTuple R counterpart of GC_lib_entry.
GC_lib_entry <- function(id, inchikey, name, RI, exact_mass, compound_formula,
                         rounded_list, peaks, base_peak, meta_text) {
  structure(list(
    id = id, inchikey = inchikey, name = name, RI = RI,
    exact_mass = exact_mass, compound_formula = compound_formula,
    rounded_list = rounded_list, peaks = as.matrix(peaks),
    base_peak = base_peak, meta_text = meta_text
  ), class = c("GC_lib_entry", "list"))
}

# Convert the two-column peak matrix into a JSON-style row list.
.gcms_matrix_rows <- function(value) {
  value <- as.matrix(value)
  if (nrow(value) == 0L) return(list())
  lapply(seq_len(nrow(value)), function(ii) as.numeric(value[ii, ]))
}

# Corresponds to load_gcms_dbfile: only accepts MSP and JSON.
load_gcms_dbfile <- function(infile) {
  if (grepl("\\.msp$", infile, ignore.case = TRUE)) {
    return(msp_standarize(parse_msp_to_listdict(infile), MSP_dict))
  }
  if (grepl("\\.json$", infile, ignore.case = TRUE)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading JSON requires jsonlite.")
    return(jsonlite::fromJSON(infile, simplifyVector = FALSE))
  }
  cat("Databse only supports MSP or JSON formats.\n")
  NULL
}

# Corresponds to reformat_gcms_lib: Convert non-uniform MSP records into fixed GC_lib_entry.
reformat_gcms_lib <- function(
    list_cpd_standards, peaks_key = "peaks", rt_key = "RETENTIONTIME",
    inchi_key = "InChIKey", name_key = "Name", mw_key = "ExactMass",
    formula_key = "Formula", filter_factor = NULL) {
  result <- list()
  counter <- 0L
  for (entry in list_cpd_standards) {
    if (!is.null(entry[[peaks_key]]) && !is.null(entry[[rt_key]]) && nzchar(entry[[rt_key]])) {
      peaks <- entry[[peaks_key]]
      base_peak <- designate_base_peak(peaks)
      if (!is.null(filter_factor) && !identical(filter_factor, FALSE)) {
        peaks <- filter_peaks_by_low_intensity_factor(peaks, base_peak[[2L]], filter_factor)
      }
      counter <- counter + 1L
      peak_matrix <- do.call(rbind, lapply(peaks, as.numeric))
      result[[length(result) + 1L]] <- GC_lib_entry(
        sprintf("cpd%08d", counter), entry[[inchi_key]] %||% "",
        entry[[name_key]] %||% "", as.numeric(entry[[rt_key]]),
        as.numeric(entry[[mw_key]]), entry[[formula_key]] %||% "",
        round(peak_matrix[, 1L]), peak_matrix, base_peak, ""
      )
    }
  }
  result
}

# The null value replacement operator makes it easy to reproduce the default value of dict.get.
`%||%` <- function(value, fallback) if (is.null(value)) fallback else value

# Corresponds to json_pseudospectra_to_msp: output MSP text record by record.
json_pseudospectra_to_msp <- function(records, outfile = "test.msp") {
  lines <- character()
  for (entry in records) {
    peaks <- entry$peaks %||% list()
    if (length(peaks) == 0L) peaks <- entry$peaks_as_features %||% list()
    lines <- c(
      lines,
      paste0("ID: ", entry$id %||% ""),
      paste0("Name: ", entry$name %||% ""),
      paste0("RETENTIONTIME: ", entry$RI %||% ""),
      paste0("InChIKey: ", entry$inchikey %||% "")
    )
    features <- entry$features %||% entry$members %||% character()
    if (length(features)) lines <- c(lines, paste0("Features: ", paste(features, collapse = ",")))
    if (nzchar(entry$annotation %||% "")) lines <- c(lines, paste0("Comment: ", entry$annotation))
    lines <- c(lines, paste0("Num Peaks: ", length(peaks)))
    for (peak in peaks) lines <- c(lines, sprintf("%.4f %.4f", peak[[1L]], peak[[2L]]))
    lines <- c(lines, "")
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to designate_base_peak: the peak where the maximum intensity value first appears.
designate_base_peak <- function(peaks) {
  peak_list <- if (is.matrix(peaks)) .gcms_matrix_rows(peaks) else peaks
  peak_list[[which.max(vapply(peak_list, `[[`, 0, 2L))]]
}

# Corresponds to three intensity filter functions, all of which retain the Python strict inequality sign.
filter_peaks_by_intensity_factor <- function(peaks, base_peak_intensity, filter_factor = 100) {
  upper <- base_peak_intensity * filter_factor
  lower <- base_peak_intensity / filter_factor
  Filter(function(peak) peak[[2L]] < upper && peak[[2L]] > lower, peaks)
}
filter_peaks_by_low_intensity_factor <- function(peaks, base_peak_intensity, filter_factor) {
  lower <- base_peak_intensity / filter_factor
  Filter(function(peak) peak[[2L]] > lower, peaks)
}
filter_features_by_low_intensity_factor <- function(features, base_peak_intensity, filter_factor) {
  lower <- base_peak_intensity / filter_factor
  Filter(function(feature) feature$peak_area > lower, features)
}

# Corresponds to filter_against_libentry: first filter candidates by rounded integer m/z.
filter_against_libentry <- function(query_spectrum, libentry) {
  selected <- query_spectrum$rounded_mzs %in% libentry$rounded_list
  list(
    query_spectrum$members[which(selected)],
    query_spectrum$peaks[selected, , drop = FALSE]
  )
}

# Corresponds to find_entries_in_rtwindow: RI distance is strictly smaller than tol.
find_entries_in_rtwindow <- function(query, gclib, tol = 30) {
  Filter(function(entry) abs(query - entry$RI) < tol, gclib)
}

# Corresponds to read_fit_KovatsIndex_rtime: Use LOWESS to fit the second column RT to the first column RI.
read_fit_KovatsIndex_rtime <- function(KovatsIndex_file, sep = "\t", frac = 0.3) {
  kovats <- utils::read.table(KovatsIndex_file, sep = sep, header = TRUE, check.names = FALSE)
  fit <- stats::lowess(x = as.numeric(kovats[[2L]]), y = as.numeric(kovats[[1L]]), f = frac)
  cbind(fit$x, fit$y)
}

# Corresponds to append_kovats_index: the asari second unit is first divided by 60 and then linearly interpolated.
append_kovats_index <- function(list_features, ri_model) {
  lapply(list_features, function(feature) {
    feature$RI <- as.numeric(stats::approx(
      ri_model[, 1L], ri_model[, 2L], xout = feature$rtime / 60, rule = 2
    )$y)
    feature
  })
}

# Corresponds to port_pseudospectrum_to_json.
port_pseudospectrum_to_json <- function(PS, normalize_intensity = FALSE) {
  peaks <- PS$peaks
  if (isTRUE(normalize_intensity)) peaks[, 2L] <- peaks[, 2L] / max(peaks[, 2L])
  list(
    id = PS$id, rtime = PS$rtime, RI = PS$RI,
    num_features = PS$num_features, members = PS$members,
    peaks = .gcms_matrix_rows(peaks), annotation = PS$annotation
  )
}

# Corresponds to serialize_annotated_empCpds: converts only two NumPy array fields.
serialize_annotated_empCpds <- function(list_empCpds) {
  lapply(list_empCpds, function(entry) {
    entry$peaks_in_lib <- .gcms_matrix_rows(entry$peaks_in_lib)
    entry$peaks_as_features <- .gcms_matrix_rows(entry$peaks_as_features)
    entry
  })
}

# Corresponds to ri_penalty_function: unilateral sigmoid RI penalty.
ri_penalty_function <- function(abs_ri_delta, min_delta = 1, max_delta = 100) {
  if (abs_ri_delta > max_delta) 0
  else if (abs_ri_delta < min_delta) 1
  else 2 / (1 + exp(exp(1) * abs_ri_delta / max_delta))
}

# Extract data.frame rows by feature id.
.gcms_feature_row <- function(feature_dataframe, id) {
  if (is.null(rownames(feature_dataframe))) stop("feature_dataframe requires feature IDs as row names.")
  as.numeric(feature_dataframe[as.character(id), , drop = TRUE])
}

# Corresponds to filter_peaks_by_penalized_distance.
filter_peaks_by_penalized_distance <- function(
    selected_features, seed_feature, feature_dataframe,
    min_ri_delta = 1, max_ri_delta = 100, feature_distance_filter = 0.5) {
  Filter(function(feature) {
    correlation <- stats::cor(
      .gcms_feature_row(feature_dataframe, feature$id),
      .gcms_feature_row(feature_dataframe, seed_feature$id),
      use = "pairwise.complete.obs"
    )
    correlation * ri_penalty_function(
      abs(feature$RI - seed_feature$RI), min_ri_delta, max_ri_delta
    ) >= feature_distance_filter
  }, selected_features)
}

# Corresponds to find_all_matches_centurion_indexed_list.
find_all_matches_centurion_indexed_list <- function(query_mz, mz_centurion_tree, limit_ppm = 5) {
  key <- as.integer(query_mz * 100)
  tolerance <- query_mz * limit_ppm * 1e-6
  result <- list()
  for (bucket in (key - 1L):(key + 1L)) {
    for (peak in mz_centurion_tree[[as.character(bucket)]] %||% list()) {
      if (abs(peak$mz - query_mz) < tolerance) result[[length(result) + 1L]] <- peak
    }
  }
  result
}

# Construct the hundredth m/z index of the GC feature.
.gcms_centurion_tree <- function(features) {
  tree <- list()
  for (feature in features) {
    key <- as.character(as.integer(feature$mz * 100))
    tree[[key]] <- c(tree[[key]], list(feature))
  }
  tree
}

# Corresponds to get_matched_features_per_cpd: returns unique feature id by absolute m/z and RI window.
get_matched_features_per_cpd <- function(cpd, mz_centurion_tree, mz_tolerance_da, ri_tolerance) {
  ids <- character()
  for (row in seq_len(nrow(cpd$peaks))) {
    peak <- cpd$peaks[row, ]
    key <- as.integer(peak[[1L]] * 100)
    for (bucket in (key - 1L):(key + 1L)) {
      candidates <- mz_centurion_tree[[as.character(bucket)]] %||% list()
      ids <- c(ids, vapply(Filter(function(feature) {
        abs(feature$mz - peak[[1L]]) < mz_tolerance_da &&
          abs(feature$RI - cpd$RI) < ri_tolerance
      }, candidates), `[[`, "", "id"))
    }
  }
  unique(ids)
}

# Corresponds to distill_correlated_features: the one with the highest correlation sum is the quant feature.
distill_correlated_features <- function(matched_feature_ids, feature_dataframe, corr_cutoff) {
  intensities <- feature_dataframe[matched_feature_ids, , drop = FALSE]
  correlations <- stats::cor(t(as.matrix(intensities)), use = "pairwise.complete.obs")
  quant <- colnames(correlations)[which.max(colSums(correlations, na.rm = TRUE))]
  selected <- matched_feature_ids[correlations[quant, matched_feature_ids] > corr_cutoff]
  list(quant, selected)
}

# Corresponds to get_seeded_pseudospectrum: builds a pseudospectrum around seed in the RI window.
get_seeded_pseudospectrum <- function(
    seed_tag, ref_ri, seed_feature, list_features, feature_dataframe,
    min_ri_delta = 1, max_ri_delta = 100, low_peak_filter_factor = 1000,
    feature_distance_filter = NULL) {
  selected <- Filter(function(feature) abs(feature$RI - ref_ri) < max_ri_delta, list_features)
  if (!is.null(low_peak_filter_factor) && !identical(low_peak_filter_factor, FALSE)) {
    selected <- filter_features_by_low_intensity_factor(
      selected, seed_feature$peak_area, low_peak_filter_factor
    )
  }
  if (!is.null(feature_distance_filter) && !identical(feature_distance_filter, FALSE)) {
    selected <- filter_peaks_by_penalized_distance(
      selected, seed_feature, feature_dataframe, min_ri_delta,
      max_ri_delta, feature_distance_filter
    )
  }
  if (length(selected) <= 1L) return(NULL)
  PseudoSpectrum(
    seed_feature$id_number, seed_feature$rtime, ref_ri,
    round(vapply(selected, `[[`, 0, "mz")), length(selected),
    vapply(selected, `[[`, "", "id"),
    cbind(
      mz = vapply(selected, `[[`, 0, "mz"),
      intensity = vapply(selected, `[[`, 0, "peak_area")
    ), seed_tag
  )
}

# Corresponds to get_spaced_top_features: select the first feature in the sorted list for each RI interval.
get_spaced_top_features <- function(list_features_sorted, ri_gap = 100) {
  selected <- list()
  for (lower in seq(1000, 4499, by = ri_gap)) {
    match <- Filter(function(feature) lower <= feature$RI && feature$RI < lower + ri_gap, list_features_sorted)
    if (length(match)) selected[[length(selected) + 1L]] <- match[[1L]]
  }
  selected
}

# Corresponds to iterative_build_pseudospectra_by_penalizeddistance.
iterative_build_pseudospectra_by_penalizeddistance <- function(
    list_features_sorted, feature_dataframe, init_core_features = character(),
    min_ri_delta = 1, max_ri_delta = 100, low_peak_filter_factor = 1000,
    feature_distance_filter = 0.5, ri_gap = 100, max_core_features = 20000) {
  core_features <- unique(as.character(init_core_features))
  spectra <- list()
  maximum <- min(max_core_features, 0.9 * length(list_features_sorted))
  while (length(core_features) <= maximum) {
    remaining <- Filter(function(feature) !(feature$id %in% core_features), list_features_sorted)
    newly_accounted <- character()
    for (seed in get_spaced_top_features(remaining, ri_gap)) {
      new <- get_seeded_pseudospectrum(
        "", seed$RI, seed, remaining, feature_dataframe,
        min_ri_delta, max_ri_delta, low_peak_filter_factor, feature_distance_filter
      )
      if (!is.null(new)) {
        spectra[[length(spectra) + 1L]] <- new
        newly_accounted <- c(newly_accounted, new$members)
      }
    }
    # If the Python source code misses the union return value, it will loop infinitely; the R version retains the intention and exits when there is no progress.
    updated <- union(core_features, newly_accounted)
    if (length(updated) == length(core_features)) break
    core_features <- updated
  }
  list(spectra, core_features)
}

# Corresponds to have_basepeak_molecularion.
have_basepeak_molecularion <- function(base_mz, mole_mz, peaks, mz_tolerance_da = 0.005) {
  peaks <- as.matrix(peaks)
  list(
    min(abs(peaks[, 1L] - base_mz)) < mz_tolerance_da,
    min(abs(peaks[, 1L] - mole_mz)) < mz_tolerance_da
  )
}

# Calls the injectable ms_entropy backend; without backend uses a JS similarity approximation that matches peak normalized intensities.
.gcms_entropy_similarity <- function(query, reference, tolerance) {
  backend <- getOption("asariR.entropy_similarity")
  if (is.function(backend)) return(as.numeric(backend(query, reference, tolerance)))
  query <- as.matrix(query)
  reference <- as.matrix(reference)
  matches <- find_matches(reference[, 1L], query[, 1L], tolerance)
  if (length(matches[[1L]]) == 0L) return(0)
  q <- query[matches[[2L]] + 1L, 2L]
  r <- reference[matches[[1L]] + 1L, 2L]
  q <- q / sum(q)
  r <- r / sum(r)
  midpoint <- (q + r) / 2
  entropy <- function(x) -sum(x[x > 0] * log(x[x > 0]))
  max(0, 1 - (entropy(midpoint) - (entropy(q) + entropy(r)) / 2) / log(2))
}

# Corresponds to batch_lib_search_score: recall by m/z and RI, refine by correlation and calculate dual scores.
batch_lib_search_score <- function(
    list_cpds, list_features, dict_features, feature_dataframe,
    mz_tolerance_da = 0.005, ri_tolerance = 30, cosine_penalty = 1,
    corr_cutoff = 0.6) {
  results <- list()
  tree <- .gcms_centurion_tree(list_features)
  for (compound in list_cpds) {
    ids <- get_matched_features_per_cpd(compound, tree, mz_tolerance_da, ri_tolerance)
    if (length(ids)) {
      distilled <- distill_correlated_features(ids, feature_dataframe, corr_cutoff)
      quant <- distilled[[1L]]
      ids <- distilled[[2L]]
      peaks <- do.call(rbind, lapply(ids, function(id) {
        c(dict_features[[id]]$mz, dict_features[[id]]$peak_area)
      }))
      entropy_score <- .gcms_entropy_similarity(peaks, compound$peaks, mz_tolerance_da)
      cosine <- cosine_similarity(
        peaks, compound$peaks, tolerance = mz_tolerance_da,
        sqrt_transform = TRUE, penalty = cosine_penalty
      )
      results[[length(results) + 1L]] <- list(
        lib_entry = compound, quant_feature = quant,
        quant_feature_RI = dict_features[[quant]]$RI,
        candidate_feature_ids = ids, entropy_score = entropy_score,
        cosine_score = cosine[[1L]], pseudo_spec = peaks,
        num_matched_peaks = cosine[[2L]]
      )
    }
  }
  cat(length(results), " matched to library compounds.\n\n", sep = "")
  results
}

# Corresponds to curate_batch_lib_search_result: sort empCpd and feature records according to any score threshold.
curate_batch_lib_search_result <- function(
    matched_results, mz_tolerance_da = 0.005,
    score_cutoff_cosine = 0.5, score_cutoff_entropy = 0.4) {
  empirical <- list()
  annotations <- list()
  counter <- 0L
  for (matched in matched_results) {
    if (matched$entropy_score >= score_cutoff_entropy || matched$cosine_score >= score_cutoff_cosine) {
      counter <- counter + 1L
      epd_id <- sprintf("empCpd_%04d", counter)
      features <- matched$candidate_feature_ids
      experimental_peaks <- matched$pseudo_spec
      library_peaks <- matched$lib_entry$peaks
      flags <- have_basepeak_molecularion(
        matched$lib_entry$base_peak[[1L]], max(library_peaks[, 1L]),
        experimental_peaks, mz_tolerance_da
      )
      empirical[[length(empirical) + 1L]] <- list(
        id = epd_id, RI = matched$quant_feature_RI,
        entropy_score = matched$entropy_score, cosine_score = matched$cosine_score,
        lib_entry_id = matched$lib_entry$id, name = matched$lib_entry$name,
        inchikey = matched$lib_entry$inchikey, features = features,
        quant_ion = matched$quant_feature, peaks_as_features = experimental_peaks,
        peaks_in_lib = library_peaks, has_basepeak = flags[[1L]],
        has_molecularion = flags[[2L]]
      )
      for (feature in features) {
        annotations[[length(annotations) + 1L]] <- list(
          feature = feature, empCpd = epd_id, quant_ion = matched$quant_feature,
          entropy_score = matched$entropy_score, cosine_score = matched$cosine_score,
          lib_entry_id = matched$lib_entry$id, name = matched$lib_entry$name,
          inchikey = matched$lib_entry$inchikey
        )
      }
    }
  }
  list(empirical, annotations)
}

# Corresponds to write_tsv_feature_anno: write 17 columns one by one annotate feature.
write_tsv_feature_anno <- function(feature_anno_list, dict_features, dict_lib_entries, outfile) {
  header <- c(
    "feature", "mz", "rtime", "RI", "empCpd", "quant_ion",
    "score_cosine", "score_entropy", "name", "inchikey", "lib_RI",
    "delta_RI", "peak_area", "cSelectivity", "peak_shape", "snr", "detection_counts"
  )
  lines <- paste(header, collapse = "\t")
  for (annotation in feature_anno_list) {
    id <- annotation$feature
    feature <- dict_features[[id]]
    lib_ri <- dict_lib_entries[[annotation$lib_entry_id]]$RI
    values <- c(
      id, feature$mz, feature$rtime, round(feature$RI, 3), annotation$empCpd,
      annotation$quant_ion, round(annotation$cosine_score, 3),
      round(annotation$entropy_score, 3), annotation$name, annotation$inchikey,
      round(lib_ri, 3), round(feature$RI - lib_ri, 3), feature$peak_area,
      feature$cSelectivity, feature$goodness_fitting, feature$snr,
      feature$detection_counts
    )
    lines <- c(lines, paste(values, collapse = "\t"))
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to write_tsv_empCpd_anno: output summary rows for each empirical compound.
write_tsv_empCpd_anno <- function(list_empCpds, dict_features, dict_lib_entries, outfile) {
  header <- c(
    "empCpd", "name", "inchikey", "formula", "lib_exact_mass",
    "score_cosine", "score_entropy", "num_matched_peaks", "quant_ion",
    "mz_quant_ion", "rtime_quant_ion(min)", "RI_quant_ion", "lib_RI",
    "delta_RI", "has_basepeak", "has_molecularion", "presence_ratio"
  )
  lines <- paste(header, collapse = "\t")
  for (epd in list_empCpds) {
    library <- dict_lib_entries[[epd$lib_entry_id]]
    feature <- dict_features[[epd$quant_ion]]
    count <- nrow(epd$peaks_as_features)
    values <- c(
      epd$id, epd$name, epd$inchikey, library$compound_formula,
      library$exact_mass, round(epd$cosine_score, 3), round(epd$entropy_score, 3),
      count, epd$quant_ion, feature$mz, round(feature$rtime / 60, 2),
      as.integer(epd$RI), as.integer(library$RI), as.integer(epd$RI - library$RI),
      epd$has_basepeak, epd$has_molecularion, round(count / nrow(library$peaks), 2)
    )
    lines <- c(lines, paste(values, collapse = "\t"))
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to get_clusters_by_ri_hcl: Ward hierarchical clustering is performed after multiplying the correlation distance by RI penalty.
get_clusters_by_ri_hcl <- function(
    feature_dataframe_sorted, dict_features,
    ri_penalty_function = get("ri_penalty_function", inherits = TRUE),
    min_delta = 1, max_delta = 100, hcl_distance_cut = 1,
    limit_fearture_num = NULL) {
  data <- feature_dataframe_sorted
  if (!is.null(limit_fearture_num) && !identical(limit_fearture_num, FALSE)) {
    data <- utils::head(data, as.integer(limit_fearture_num))
  }
  number <- nrow(data)
  cat("Working with ", number, " features.\n", sep = "")
  if (number <= 1L) {
    clusters <- rep(1L, number)
  } else {
    ids <- rownames(data)
    ri <- vapply(ids, function(id) dict_features[[id]]$RI, 0)
    correlation <- stats::cor(t(as.matrix(data)), use = "pairwise.complete.obs")
    distance <- matrix(0, number, number)
    for (ii in seq_len(number - 1L)) {
      for (jj in (ii + 1L):number) {
        penalty <- ri_penalty_function(abs(ri[[ii]] - ri[[jj]]), min_delta, max_delta)
        value <- 1 - penalty * correlation[ii, jj]
        distance[ii, jj] <- distance[jj, ii] <- value
      }
    }
    tree <- stats::hclust(stats::as.dist(distance), method = "ward.D2")
    clusters <- stats::cutree(tree, h = hcl_distance_cut)
  }
  cat("number of features: ", length(clusters), "\n")
  cat("number of clusters: ", length(unique(clusters)), "\n")
  cluster_dict <- split(rownames(data), clusters)
  list(as.integer(clusters), unname(cluster_dict))
}

# Corresponds to extend_cluster: use the highest area feature of the original cluster as the seed to extend the relevant members.
extend_cluster <- function(
    fcluster, featureDict, list_features, feature_dataframe,
    ri_tolerance = 50, correlation_cut = 0.7) {
  cluster_features <- lapply(fcluster, function(id) featureDict[[id]])
  areas <- vapply(cluster_features, `[[`, 0, "peak_area")
  seed <- cluster_features[[order(-areas)[[1L]]]]
  candidates <- Filter(function(feature) abs(feature$RI - seed$RI) < ri_tolerance, list_features)
  selected <- vapply(Filter(function(feature) {
    stats::cor(
      .gcms_feature_row(feature_dataframe, feature$id),
      .gcms_feature_row(feature_dataframe, seed$id),
      use = "pairwise.complete.obs"
    ) >= correlation_cut
  }, candidates), `[[`, "", "id")
  list(c(fcluster, selected), seed)
}

# Corresponds to format_fcluster_to_pseudospectrum.
format_fcluster_to_pseudospectrum <- function(selected_features, seed_feature) {
  PseudoSpectrum(
    seed_feature$id_number, seed_feature$rtime, seed_feature$RI,
    round(vapply(selected_features, `[[`, 0, "mz")), length(selected_features),
    vapply(selected_features, `[[`, "", "id"),
    cbind(
      mz = vapply(selected_features, `[[`, 0, "mz"),
      intensity = vapply(selected_features, `[[`, 0, "peak_area")
    ), ""
  )
}

# Corresponds to iterative_build_pseudospectra_by_hcl.
iterative_build_pseudospectra_by_hcl <- function(
    list_features_sorted, feature_dataframe, init_core_features = character(),
    step_size = 2000, hcl_distance_cut = 1, ri_tolerance = 50,
    correlation_cut = 0.7, max_core_features = 20000) {
  core <- unique(as.character(init_core_features))
  spectra <- list()
  maximum <- min(max_core_features, 0.9 * length(list_features_sorted))
  feature_dict <- stats::setNames(list_features_sorted, vapply(list_features_sorted, `[[`, "", "id"))
  counter <- 0L
  while (length(core) <= maximum) {
    counter <- counter + 1L
    remaining <- Filter(function(feature) !(feature$id %in% core), list_features_sorted)
    cat("Iteration ", counter, " of building pseudospectra from ", length(remaining), " remaining features.\n", sep = "")
    if (length(remaining) < 2L) break
    clustered <- get_clusters_by_ri_hcl(
      feature_dataframe[vapply(remaining, `[[`, "", "id"), , drop = FALSE],
      feature_dict, hcl_distance_cut = hcl_distance_cut,
      limit_fearture_num = step_size
    )[[2L]]
    tail_features <- if (length(remaining) > step_size) remaining[-seq_len(step_size)] else list()
    new_clusters <- lapply(clustered[lengths(clustered) > 1L], function(cluster) {
      extend_cluster(cluster, feature_dict, tail_features, feature_dataframe, ri_tolerance, correlation_cut)
    })
    accounted <- character()
    for (cluster in new_clusters) {
      spectra[[length(spectra) + 1L]] <- format_fcluster_to_pseudospectrum(
        lapply(cluster[[1L]], function(id) feature_dict[[id]]), cluster[[2L]]
      )
      accounted <- c(accounted, cluster[[1L]])
    }
    updated <- union(core, accounted)
    if (length(updated) == length(core)) break
    core <- updated
    cat("Core features collected: ", length(core), ".\n", sep = "")
  }
  list(spectra, core)
}

# Corresponds to batch_lib_search_by_basepeaks: Use the library base peak as the seed to create and score pseudo-spectra.
batch_lib_search_by_basepeaks <- function(
    list_lib_entries, list_features, feature_dataframe,
    ms2_tolerance_in_ppm = 5, ms2_tolerance_in_da = 0.005,
    ri_tolerance = 50, cosine_penalty = 1, min_ri_delta = 1,
    max_ri_delta = 100, low_peak_filter_factor = 1000,
    feature_distance_filter = NULL) {
  results <- list()
  compound_dict <- stats::setNames(list_lib_entries, vapply(list_lib_entries, `[[`, "", "id"))
  feature_dict <- stats::setNames(list_features, vapply(list_features, `[[`, "", "id"))
  compounds <- lapply(list_lib_entries, function(entry) {
    list(id = entry$id, mz = entry$base_peak[[1L]], rtime = entry$RI)
  })
  features <- lapply(list_features, function(feature) {
    list(id = feature$id, mz = feature$mz, rtime = feature$RI)
  })
  matched <- list_match_lcms_features(
    compounds, features, ms2_tolerance_in_ppm, ri_tolerance
  )
  for (library_id in names(matched)) {
    for (feature_id in matched[[library_id]]) {
      entry <- compound_dict[[library_id]]
      pseudo <- get_seeded_pseudospectrum(
        library_id, entry$RI, feature_dict[[feature_id]], list_features,
        feature_dataframe, min_ri_delta, max_ri_delta,
        low_peak_filter_factor, feature_distance_filter
      )
      if (!is.null(pseudo)) {
        filtered <- filter_against_libentry(pseudo, entry)
        entropy_score <- .gcms_entropy_similarity(filtered[[2L]], entry$peaks, ms2_tolerance_in_da)
        cosine <- cosine_similarity(
          filtered[[2L]], entry$peaks, tolerance = ms2_tolerance_in_da,
          sqrt_transform = TRUE, penalty = cosine_penalty
        )
        results[[length(results) + 1L]] <- list(
          lib_entry = entry, pseudo_spec = pseudo,
          candidate_features = filtered[[1L]], base_feature_id = feature_id,
          entropy_score = entropy_score, cosine_score = cosine[[1L]],
          num_matched_peaks = cosine[[2L]]
        )
      }
    }
  }
  cat(length(results), " matched results found by base peak search.\n", sep = "")
  results
}

# Corresponds to export_feature_annotation_bybasepeaksearch.
export_feature_annotation_bybasepeaksearch <- function(
    matched_results, feature_dataframe, score_cutoff_cosine = 0.5,
    score_cutoff_entropy = 0.4, corr_cutoff = 0.7, mz_tolerance_ppm = 5) {
  empirical <- list()
  annotations <- list()
  counter <- 0L
  for (matched in matched_results) {
    mapping <- all_mass_paired_mapping(
      matched$lib_entry$peaks[, 1L], matched$pseudo_spec$peaks[, 1L],
      std_ppm = mz_tolerance_ppm
    )$mapped
    if (length(mapping)) {
      counter <- counter + 1L
      epd_id <- sprintf("empCpd_%04d", counter)
      experimental_idx <- vapply(mapping, `[[`, 0L, 2L)
      library_idx <- vapply(mapping, `[[`, 0L, 1L)
      features <- matched$pseudo_spec$members[experimental_idx]
      experimental_peaks <- matched$pseudo_spec$peaks[experimental_idx, , drop = FALSE]
      library_peaks <- matched$lib_entry$peaks[library_idx, , drop = FALSE]
      quant <- matched$base_feature_id
      if (matched$entropy_score >= score_cutoff_entropy || matched$cosine_score >= score_cutoff_cosine) {
        empirical[[length(empirical) + 1L]] <- list(
          id = epd_id, RI = matched$pseudo_spec$RI,
          entropy_score = matched$entropy_score, cosine_score = matched$cosine_score,
          lib_entry_id = matched$lib_entry$id, name = matched$lib_entry$name,
          inchikey = matched$lib_entry$inchikey, features = features,
          quant_ion = quant, peaks_as_features = experimental_peaks,
          peaks_in_lib = library_peaks
        )
        for (feature in features) {
          correlation <- stats::cor(
            .gcms_feature_row(feature_dataframe, quant),
            .gcms_feature_row(feature_dataframe, feature),
            use = "pairwise.complete.obs"
          )
          annotations[[length(annotations) + 1L]] <- list(
            feature = feature, empCpd = epd_id, quant_ion = quant,
            entropy_score = matched$entropy_score, cosine_score = matched$cosine_score,
            lib_entry_id = matched$lib_entry$id, name = matched$lib_entry$name,
            inchikey = matched$lib_entry$inchikey, correlation = correlation,
            is_core = correlation >= corr_cutoff
          )
        }
      }
    }
  }
  list(empirical, annotations)
}

# Corresponds to group_pseudospectra_from_features: Use the highest area feature to create RT window pseudospectra in sequence.
group_pseudospectra_from_features <- function(
    list_features, rtime_window_in_seconds = 1, bin_fraction = 0.2) {
  spectra <- list()
  counted <- character()
  number_seed <- as.integer(length(list_features) * bin_fraction)
  if (number_seed > 0L) {
    for (feature in list_features[seq_len(number_seed)]) {
      if (!(as.character(feature$id_number) %in% counted)) {
        in_range <- Filter(function(candidate) {
          abs(candidate$rtime - feature$rtime) < rtime_window_in_seconds
        }, list_features)
        # Keep the original Python implementation: peaks reuses mz and peak_area of seed for each row.
        peaks <- matrix(
          rep(c(feature$mz, feature$peak_area), length(in_range)),
          ncol = 2L, byrow = TRUE
        )
        spectra[[length(spectra) + 1L]] <- PseudoSpectrum(
          feature$id_number, feature$rtime, feature$RI,
          rep(round(feature$mz), length(in_range)), length(in_range),
          vapply(in_range, function(item) as.character(item$id_number), ""),
          peaks, ""
        )
        counted <- union(counted, vapply(in_range, function(item) as.character(item$id_number), ""))
      }
    }
  }
  cat(sprintf(
    "From %d features, %d pseudo spectra were constructed, of which %d have 100 or more peaks.\n",
    length(list_features), length(spectra),
    sum(vapply(spectra, function(item) item$num_features >= 100L, FALSE))
  ))
  spectra
}

# Corresponds to reverse_spec_searches: filter by entropy and reverse cosine respectively in the RI window.
reverse_spec_searches <- function(
    list_pseudo_spectra, list_lib_entries, ri_window = 100,
    mz_tolerance = 0.005, cosine_penalty = 1, score_cutoff = 0.5,
    score_cutoff_entropy = 0.4) {
  matched_cosine <- list()
  matched_entropy <- list()
  for (spectrum in list_pseudo_spectra) {
    candidates <- find_entries_in_rtwindow(spectrum$RI, list_lib_entries, ri_window)
    for (entry in candidates) {
      filtered <- filter_against_libentry(spectrum, entry)
      entropy_score <- .gcms_entropy_similarity(filtered[[2L]], entry$peaks, mz_tolerance)
      if (entropy_score > score_cutoff_entropy) {
        matched_entropy[[length(matched_entropy) + 1L]] <- list(
          entropy_score, entry, filtered[[1L]], spectrum
        )
      }
      cosine <- cosine_similarity(
        filtered[[2L]], entry$peaks, tolerance = mz_tolerance,
        sqrt_transform = TRUE, penalty = cosine_penalty
      )
      if (cosine[[1L]] > score_cutoff) {
        matched_cosine[[length(matched_cosine) + 1L]] <- list(
          cosine[[1L]], entry, filtered[[1L]], spectrum, cosine[[2L]]
        )
      }
    }
  }
  cat(sprintf(
    "Found %d by entropy serach and %d by cosine search. \n",
    length(matched_entropy), length(matched_cosine)
  ))
  list(matched_entropy, matched_cosine)
}

# Corresponds to export_feature_annotations: output simple feature matching tuple.
export_feature_annotations <- function(matched_list, mz_tolerance_ppm = 5) {
  annotations <- list()
  for (matched in matched_list) {
    mapping <- complete_mass_paired_mapping(
      matched[[2L]]$peaks[, 1L], matched[[4L]]$peaks[, 1L],
      std_ppm = mz_tolerance_ppm
    )$mapped
    for (pair in mapping) {
      annotations[[length(annotations) + 1L]] <- list(
        matched[[4L]]$members[[pair[[2L]]]], matched[[1L]],
        matched[[2L]]$name, matched[[2L]]$inchikey, matched[[4L]]$RI
      )
    }
  }
  annotations
}

# Corresponds to export_feature_annotation_details: select the quant feature with the highest joint strength ranking.
export_feature_annotation_details <- function(
    matched_list, feature_dataframe, corr_cutoff = 0.7, mz_tolerance_ppm = 5) {
  empirical <- list()
  annotations <- list()
  counter <- 0L
  for (matched in matched_list) {
    mapping <- complete_mass_paired_mapping(
      matched[[2L]]$peaks[, 1L], matched[[4L]]$peaks[, 1L],
      std_ppm = mz_tolerance_ppm
    )$mapped
    if (length(mapping)) {
      counter <- counter + 1L
      epd_id <- sprintf("empCpd_%04d", counter)
      exp_idx <- vapply(mapping, `[[`, 0L, 2L)
      lib_idx <- vapply(mapping, `[[`, 0L, 1L)
      features <- matched[[4L]]$members[exp_idx]
      exp_peaks <- matched[[4L]]$peaks[exp_idx, , drop = FALSE]
      lib_peaks <- matched[[2L]]$peaks[lib_idx, , drop = FALSE]
      top <- which.max((order(exp_peaks[, 2L]) - 1L) + (order(lib_peaks[, 2L]) - 1L))
      quant <- features[[top]]
      empirical[[length(empirical) + 1L]] <- list(
        id = epd_id, RI = matched[[4L]]$RI, score = matched[[1L]],
        name = matched[[2L]]$name, inchikey = matched[[2L]]$inchikey,
        features = features, quant_ion = quant,
        peaks_as_features = exp_peaks, peaks_in_lib = lib_peaks
      )
      for (feature in features) {
        correlation <- stats::cor(
          .gcms_feature_row(feature_dataframe, quant),
          .gcms_feature_row(feature_dataframe, feature),
          use = "pairwise.complete.obs"
        )
        annotations[[length(annotations) + 1L]] <- list(
          feature = feature, empCpd = epd_id, quant_ion = quant,
          score = matched[[1L]], name = matched[[2L]]$name,
          inchikey = matched[[2L]]$inchikey, correlation = correlation,
          is_core = correlation >= corr_cutoff
        )
      }
    }
  }
  list(empirical, annotations)
}

# Corresponds to iterative_reverse_annotation: multiple rounds of excluding core features and merging the two scoring results.
iterative_reverse_annotation <- function(
    list_features, list_lib_entries, feature_dataframe,
    binning_rtime_window_in_seconds = 1, search_ri_window = 50,
    search_mz_tolerance = 0.005, cosine_penalty = 1,
    score_cutoff = 0.5, score_cutoff_entropy = 0.4, corr_cutoff = 0.7,
    export_mz_tolerance_ppm = 5, bin_fraction = 0.2, iterations = 3L) {
  core <- character()
  epd_cosine <- anno_cosine <- epd_entropy <- anno_entropy <- list()
  for (step in seq_len(iterations)) {
    cat("Iteration step ", step, "\n")
    fraction <- if (step == iterations) 1 else bin_fraction
    spectra <- group_pseudospectra_from_features(
      Filter(function(feature) !(feature$id %in% core), list_features),
      binning_rtime_window_in_seconds, fraction
    )
    matched <- reverse_spec_searches(
      spectra, list_lib_entries, search_ri_window, search_mz_tolerance,
      cosine_penalty, score_cutoff, score_cutoff_entropy
    )
    cosine <- export_feature_annotation_details(
      matched[[2L]], feature_dataframe, corr_cutoff, export_mz_tolerance_ppm
    )
    entropy <- export_feature_annotation_details(
      matched[[1L]], feature_dataframe, corr_cutoff, export_mz_tolerance_ppm
    )
    epd_cosine <- c(epd_cosine, cosine[[1L]])
    anno_cosine <- c(anno_cosine, cosine[[2L]])
    epd_entropy <- c(epd_entropy, entropy[[1L]])
    anno_entropy <- c(anno_entropy, entropy[[2L]])
    core <- union(core, vapply(
      Filter(function(item) isTRUE(item$is_core), c(anno_cosine, anno_entropy)),
      `[[`, "", "feature"
    ))
  }
  list(
    list_empCpds_cosine = epd_cosine,
    feature_anno_list_cosine = anno_cosine,
    list_empCpds_entropy = epd_entropy,
    feature_anno_list_entropy = anno_entropy,
    core_features = core
  )
}

# Corresponds to two equality functions.
is_same_match <- function(epd1, epd2) {
  identical(epd1$inchikey, epd2$inchikey) &&
    identical(epd1$quant_ion, epd2$quant_ion) &&
    identical(epd1$features, epd2$features)
}
is_same_matched_feature <- function(f1, f2) {
  identical(f1$inchikey, f2$inchikey) && identical(f1$quant_ion, f2$quant_ion)
}

# Corresponds to the merge function: update the old record with the new record field when matching, otherwise append.
merge_epd_with_list <- function(epd1, list2_epds) {
  matched <- FALSE
  for (ii in seq_along(list2_epds)) {
    if (is_same_match(list2_epds[[ii]], epd1)) {
      list2_epds[[ii]][names(epd1)] <- epd1
      matched <- TRUE
    }
  }
  if (!matched) list2_epds[[length(list2_epds) + 1L]] <- epd1
  list2_epds
}
merge_feature_with_list <- function(feat, list2_features) {
  matched <- FALSE
  for (ii in seq_along(list2_features)) {
    if (is_same_matched_feature(feat, list2_features[[ii]])) {
      list2_features[[ii]][names(feat)] <- feat
      matched <- TRUE
    }
  }
  if (!matched) list2_features[[length(list2_features) + 1L]] <- feat
  list2_features
}

# Corresponds to format_arrays_to_lists.
format_arrays_to_lists <- function(epd) {
  epd$peaks_as_features <- .gcms_matrix_rows(epd$peaks_as_features)
  epd$peaks_in_lib <- .gcms_matrix_rows(epd$peaks_in_lib)
  epd
}

# Corresponds to cleanup_anno_empcpds_features: merge duplicate records of cosine and entropy algorithms.
cleanup_anno_empcpds_features <- function(resultDict) {
  empirical <- list()
  feature_annotations <- list()
  for (epd in resultDict$list_empCpds_cosine) {
    epd <- format_arrays_to_lists(epd)
    key <- epd$inchikey
    epd$score_cosine <- epd$score
    epd$score <- NULL
    existing <- empirical[[key]] %||% list()
    if (!any(vapply(existing, is_same_match, FALSE, epd2 = epd))) existing[[length(existing) + 1L]] <- epd
    empirical[[key]] <- existing
  }
  for (epd in resultDict$list_empCpds_entropy) {
    epd <- format_arrays_to_lists(epd)
    key <- epd$inchikey
    epd$score_entropy <- epd$score
    epd$score <- NULL
    empirical[[key]] <- if (is.null(empirical[[key]])) list(epd) else merge_epd_with_list(epd, empirical[[key]])
  }
  for (feature in resultDict$feature_anno_list_cosine) {
    key <- feature$feature
    feature$score_cosine <- feature$score
    feature$score <- NULL
    existing <- feature_annotations[[key]] %||% list()
    if (!any(vapply(existing, is_same_matched_feature, FALSE, f2 = feature))) existing[[length(existing) + 1L]] <- feature
    feature_annotations[[key]] <- existing
  }
  for (feature in resultDict$feature_anno_list_entropy) {
    key <- feature$feature
    feature$score_entropy <- feature$score
    feature$score <- NULL
    feature_annotations[[key]] <- if (is.null(feature_annotations[[key]])) {
      list(feature)
    } else merge_feature_with_list(feature, feature_annotations[[key]])
  }
  list(empirical, feature_annotations)
}

# Corresponds to append_anno_json: collected by feature and sorted in reverse (score, name, inchikey).
append_anno_json <- function(list_features, search_matched_features) {
  grouped <- list()
  for (entry in search_matched_features) {
    key <- as.character(entry[[1L]])
    grouped[[key]] <- c(grouped[[key]], list(list(round(entry[[2L]], 3), entry[[3L]], entry[[4L]])))
  }
  lapply(list_features, function(feature) {
    matches <- grouped[[as.character(feature$id)]]
    if (length(matches)) {
      keys <- vapply(matches, function(x) paste(sprintf("%020.3f", x[[1L]]), x[[2L]], x[[3L]], sep = "\r"), "")
      feature$annotation <- matches[order(keys, decreasing = TRUE)]
    }
    feature
  })
}

# Corresponds to export_jsonanno_tsv: after adding top_match/all_matches, the output is fixed at 14 columns.
export_jsonanno_tsv <- function(list_features, outfile) {
  header <- c(
    "id_number", "mz", "RI", "top_match", "all_matches", "rtime",
    "left_base", "right_base", "parent_masstrack_id", "peak_area",
    "cSelectivity", "goodness_fitting", "snr", "detection_counts"
  )
  lines <- paste(header, collapse = "\t")
  for (feature in list_features) {
    feature$top_match <- ""
    feature$all_matches <- list()
    if (!is.null(feature$annotation)) {
      feature$top_match <- feature$annotation[[1L]][[2L]]
      feature$all_matches <- feature$annotation
    }
    values <- vapply(header, function(field) {
      value <- feature[[field]]
      if (is.list(value)) paste(capture.output(dput(value)), collapse = "") else as.character(value)
    }, "")
    lines <- c(lines, paste(values, collapse = "\t"))
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}
