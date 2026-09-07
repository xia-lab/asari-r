# Corresponds to Python asari/tools/ms2.py: extract, filter and match MS2 spectra to MS1 features.

# Unify the test injected spectra or mzML files into a standard scan list.
.ms2_read_spectra <- function(infile) {
  if (is.list(infile)) return(infile)
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Reading mzML requires the optional Bioconductor package mzR.", call. = FALSE)
  }
  handle <- mzR::openMSfile(infile)
  on.exit(mzR::close(handle), add = TRUE)
  header <- mzR::header(handle)
  lapply(seq_len(nrow(header)), function(ii) {
    precursor <- if ("precursorMZ" %in% names(header)) header$precursorMZ[[ii]] else NA_real_
    list(
      ms_level = as.integer(header$msLevel[[ii]]),
      precursor_mz = if (is.na(precursor)) NULL else as.numeric(precursor),
      rtime = as.numeric(header$retentionTime[[ii]]),
      peaks = mzR::peaks(handle, ii)
    )
  })
}

# Normalize a matrix or peak tuple list into a two-column numerical matrix.
.ms2_peak_matrix <- function(peaks) {
  if (is.null(peaks) || length(peaks) == 0L) return(matrix(numeric(), ncol = 2L))
  result <- if (is.matrix(peaks)) peaks[, 1:2, drop = FALSE] else do.call(rbind, lapply(peaks, as.numeric))
  storage.mode(result) <- "double"
  result
}

# Corresponds to extract_all_spectra_form_file: extract MS1, MS2 and other level spectra at one time.
extract_all_spectra_form_file <- function(infile, min_intensity = 1000, MS2_peak_limit = 50L) {
  ms1 <- list()
  ms2 <- list()
  others <- list()
  last_precursor <- NULL
  spectra <- .ms2_read_spectra(infile)
  for (counter in seq_along(spectra)) {
    spectrum <- spectra[[counter]]
    peaks <- .ms2_peak_matrix(spectrum$peaks)
    level <- as.integer(spectrum$ms_level)
    if (level == 2L) {
      precursor <- spectrum$precursor_mz
      last_precursor <- precursor
      if (is.null(precursor) || is.na(precursor)) next
      keep <- peaks[, 2L] >= min_intensity & peaks[, 1L] < precursor - 1
      selected <- peaks[keep, , drop = FALSE]
      if (nrow(selected) > 0L) {
        selected[, 2L] <- as.integer(selected[, 2L])
        if (!is.null(MS2_peak_limit) && !identical(MS2_peak_limit, FALSE)) {
          selected <- selected[order(-selected[, 2L]), , drop = FALSE]
          selected <- utils::head(selected, as.integer(MS2_peak_limit))
        }
        ms2[[length(ms2) + 1L]] <- list(
          id = paste0("sp", counter),
          precursor_mz = precursor,
          rtime = spectrum$rtime,
          peaks = lapply(seq_len(nrow(selected)), function(ii) as.numeric(selected[ii, ]))
        )
      }
    } else if (level == 1L) {
      selected <- peaks[peaks[, 2L] >= min_intensity, , drop = FALSE]
      if (nrow(selected) > 0L) selected[, 2L] <- as.integer(selected[, 2L])
      ms1[[length(ms1) + 1L]] <- list(
        id = paste0("sp", counter),
        rtime = spectrum$rtime,
        peaks = lapply(seq_len(nrow(selected)), function(ii) as.numeric(selected[ii, ]))
      )
    } else {
      selected <- peaks[peaks[, 2L] >= min_intensity, , drop = FALSE]
      if (nrow(selected) > 0L) selected[, 2L] <- as.integer(selected[, 2L])
      others[[length(others) + 1L]] <- list(
        ms_level = level,
        id = paste0("sp", counter),
        precursor_mz = last_precursor,
        rtime = spectrum$rtime,
        peaks = lapply(seq_len(nrow(selected)), function(ii) as.numeric(selected[ii, ]))
      )
    }
  }
  list(ms1, ms2, others)
}

# Corresponds to regroup_matches_per_file: concentrate spectral matching by database INCHIKEY.
regroup_matches_per_file <- function(results) {
  grouped <- list()
  tally <- list()
  for (record in results) {
    key <- record[[1L]]$INCHIKEY
    tally[[key]] <- record[[1L]]
    remainder <- record[-1L]
    if (!is.null(grouped[[key]])) {
      grouped[[key]]$matched_spectra[[length(grouped[[key]]$matched_spectra) + 1L]] <- remainder
    } else {
      grouped[[key]] <- list(
        DB_record = list(INCHIKEY = key, name = record[[1L]]$NAME),
        matched_spectra = list(remainder)
      )
    }
  }
  list(grouped, tally)
}

# Corresponds to get_ms1_match_to_ms2cluster: select the highest scoring MS2 and the strongest MS1 peak in the window.
get_ms1_match_to_ms2cluster <- function(db_matched_cluster, list_ms1_spectra, mz_tol_ppm = 5, rt_tol = 1) {
  matches <- db_matched_cluster$matched_spectra
  scores <- vapply(matches, `[[`, 0, 4L)
  best_ms2 <- matches[[order(-scores)[[1L]]]]
  spectra <- Filter(function(spectrum) abs(spectrum$rtime - best_ms2[[3L]]) <= rt_tol, list_ms1_spectra)
  mz_tol <- mz_tol_ppm * 1e-6 * best_ms2[[2L]]
  peaks <- unlist(lapply(spectra, function(spectrum) {
    Filter(function(peak) abs(peak[[1L]] - best_ms2[[2L]]) < mz_tol, spectrum$peaks)
  }), recursive = FALSE)
  if (length(peaks) == 0L) return(list(best_ms2, NULL))
  intensities <- vapply(peaks, `[[`, 0, 2L)
  list(best_ms2, peaks[[order(intensities)[[length(intensities)]]]])
}

# Corresponds to export_table_ms1match_results: output the best MS1 peak height for each compound and each sample.
export_table_ms1match_results <- function(resultDict, cpdDict, outfile) {
  samples <- names(resultDict)
  header1 <- c("NAME", "INCHIKEY", "precursor_mz", "PRECURSORTYPE", "SMILES", "Num Peaks")
  header2 <- sub("_MS2.json", "", samples, fixed = TRUE)
  lines <- paste(c(header1, header2), collapse = "\t")
  for (key in names(cpdDict)) {
    compound <- cpdDict[[key]]
    part1 <- vapply(header1, function(field) as.character(compound[[field]]), "")
    part2 <- rep("", length(samples))
    for (ii in seq_along(samples)) {
      in_sample <- resultDict[[samples[[ii]]]][[key]]
      if (!is.null(in_sample) && !is.null(in_sample$best_ms1_peak)) {
        part2[[ii]] <- as.character(in_sample$best_ms1_peak[[2L]])
      }
    }
    lines <- c(lines, paste(c(part1, part2), collapse = "\t"))
  }
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to extract_ms2_form_file: extract only MS2 scans and remove peaks near precursor ions.
extract_ms2_form_file <- function(infile, min_intensity = 1000) {
  result <- list()
  for (spectrum in .ms2_read_spectra(infile)) {
    if (as.integer(spectrum$ms_level) != 2L) next
    precursor <- spectrum$precursor_mz
    if (is.null(precursor) || is.na(precursor)) next
    peaks <- .ms2_peak_matrix(spectrum$peaks)
    selected <- peaks[peaks[, 2L] >= min_intensity & peaks[, 1L] < precursor - 1, , drop = FALSE]
    if (nrow(selected) > 0L) {
      result[[length(result) + 1L]] <- list(
        precursor_mz = precursor,
        rtime = spectrum$rtime,
        peaks = lapply(seq_len(nrow(selected)), function(ii) as.numeric(selected[ii, ]))
      )
    }
  }
  result
}

# Corresponds to get_top_n_ms2_spectra: select the top n spectra in descending order of peak number.
get_top_n_ms2_spectra <- function(ms2_spectra, n = 5L) {
  if (length(ms2_spectra) <= n) return(ms2_spectra)
  counts <- vapply(ms2_spectra, function(spectrum) length(spectrum$peaks), 0L)
  ms2_spectra[order(-counts, seq_along(counts))[seq_len(n)]]
}

# Corresponds to get_top_n_peaks: select the top n peaks in descending order of intensity.
get_top_n_peaks <- function(spectrum, n = 50L) {
  if (length(spectrum) <= n) return(spectrum)
  intensity <- vapply(spectrum, `[[`, 0, 2L)
  spectrum[order(-intensity, seq_along(intensity))[seq_len(n)]]
}

# Corresponds to get_best_ms2_spectrum: select the best spectrum based on the total intensity of all peaks.
get_best_ms2_spectrum <- function(ms2_spectra) {
  if (length(ms2_spectra) == 0L) return(NULL)
  totals <- vapply(ms2_spectra, function(spectrum) {
    sum(vapply(spectrum$peaks, `[[`, 0, 2L))
  }, 0)
  ms2_spectra[[order(-totals, seq_along(totals))[[1L]]]]
}

# Corresponds to get_matched_ms2_ms1: first match by m/z, then filter by RT and select the best MS2 for each MS1.
get_matched_ms2_ms1 <- function(LCMS1_features, list_ms2_spectra, rt_tol = 30, ppm_tol = 5) {
  mapping <- complete_mass_paired_mapping(
    vapply(LCMS1_features, `[[`, 0, "mz"),
    vapply(list_ms2_spectra, `[[`, 0, "precursor_mz"),
    std_ppm = ppm_tol
  )$mapped
  mapping <- Filter(function(pair) {
    abs(LCMS1_features[[pair[[1L]]]]$rtime - list_ms2_spectra[[pair[[2L]]]]$rtime) <= rt_tol
  }, mapping)
  grouped <- list()
  for (pair in mapping) {
    key <- as.character(pair[[1L]])
    grouped[[key]] <- c(grouped[[key]], list(list_ms2_spectra[[pair[[2L]]]]))
  }
  lapply(names(grouped), function(key) {
    spectra <- grouped[[key]]
    list(LCMS1_features[[as.integer(key)]]$id, list(get_best_ms2_spectrum(spectra), length(spectra)))
  })
}

# Corresponds to match_ms2files_to_features: accumulate the best MS2 file by file and write out JSON.
match_ms2files_to_features <- function(
    ms1_fulltable,
    list_ms2_files,
    rt_tol = 30,
    ppm_tol = 5,
    output_json = "matched_ms2_spectra.json") {
  LCMS_features <- read_features_from_asari_table(
    paste(readLines(ms1_fulltable, warn = FALSE), collapse = "\n")
  )[[2L]]
  master <- list()
  for (ms2_file in list_ms2_files) {
    cat("Processing MS2 file: ", as.character(ms2_file), "\n", sep = "")
    spectra <- extract_ms2_form_file(ms2_file)
    cat("  Extracted ", length(spectra), " MS2 spectra\n", sep = "")
    matched <- get_matched_ms2_ms1(LCMS_features, spectra, rt_tol, ppm_tol)
    cat("  Found ", length(matched), " features with matched MS2 spectra\n", sep = "")
    for (item in matched) {
      feature_id <- item[[1L]]
      spectrum <- item[[2L]][[1L]]
      count <- item[[2L]][[2L]]
      if (is.null(master[[feature_id]])) {
        master[[feature_id]] <- list(spectrum, count)
      } else {
        total <- master[[feature_id]][[2L]] + count
        best <- get_best_ms2_spectrum(list(master[[feature_id]][[1L]], spectrum))
        master[[feature_id]] <- list(best, total)
      }
    }
  }
  # Python's current implementation passes in the entire spectrum dict and leaves it unchanged if the length is less than 20.
  for (feature_id in names(master)) master[[feature_id]] <- list(master[[feature_id]][[1L]], master[[feature_id]][[2L]])
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Writing MS2 JSON requires the optional jsonlite package.", call. = FALSE)
  }
  jsonlite::write_json(master, output_json, pretty = TRUE, auto_unbox = TRUE)
  master
}
