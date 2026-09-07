# Corresponds to Python asari/tools/qc.py: feature table QC diagram and single mzML file HTML report.

# The NULL backup value of this module ensures that it can also run with a separate source.
.qc_or <- function(value, fallback) if (is.null(value)) fallback else value

# Corresponds to asari_qc_plot: log2 area/SNR, peak shape, cSelectivity and point size combination plot.
asari_qc_plot <- function(data, outfile = "qc_plot.pdf", height = 12,
                          aspect = 0.7, cmap = NULL) {
  invisible(cmap)
  data$log2snr <- log2(data$snr)
  data$log2area <- log2(data$peak_area)
  grDevices::pdf(outfile, width = height * aspect, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  colors <- grDevices::hcl.colors(100, "Spectral")[pmax(1L, pmin(100L, ceiling(data$cSelectivity * 100)))]
  sizes <- 0.3 + 2.7 * (data$log2snr - min(data$log2snr)) /
    max(1e-12, max(data$log2snr) - min(data$log2snr))
  graphics::plot(
    data$log2area, data$goodness_fitting, pch = 16, col = colors, cex = sizes,
    ylim = c(0, 1.05), xlab = "log2area", ylab = "goodness_fitting"
  )
  graphics::grid()
  invisible(NULL)
}

# Corresponds to generate_qc_report internal calcTIC.
.qc_calcTIC <- function(spectra, mslevel = 1L) {
  rt <- numeric(); tic <- numeric(); scan_no <- 0L
  for (spectrum in spectra) if (spectrum$ms_level == mslevel) {
    scan_no <- scan_no + 1L
    peaks <- .ms2_peak_matrix(spectrum$peaks)
    tic <- c(tic, sum(as.integer(peaks[, 2L])))
    rt <- c(rt, scan_no)
  }
  list(rt, tic)
}

# Corresponds to generate_qc_report internal extract_trio.
.qc_extract_trio <- function(spectra, mslevel = 1L, min_intensity = 1000) {
  result <- list()
  for (spectrum in spectra) if (spectrum$ms_level == mslevel) {
    peaks <- .ms2_peak_matrix(spectrum$peaks)
    for (ii in which(as.integer(peaks[, 2L]) > min_intensity)) {
      result[[length(result) + 1L]] <- c(as.integer(peaks[ii, 2L]), peaks[ii, 1L], spectrum$rtime)
    }
  }
  result
}

# Corresponds to find_targets inside generate_qc_report.
.qc_find_targets <- function(trio_list, spikeins, mz_error = 0.003, min_intensity = 10000) {
  result <- list()
  for (target in spikeins) {
    result[[as.character(target[[1L]])]] <- Filter(function(trio) {
      trio[[1L]] > min_intensity && target[[2L]] - mz_error < trio[[2L]] && trio[[2L]] < target[[2L]] + mz_error
    }, trio_list)
  }
  result
}

# Simple HTML escaping.
.qc_html <- function(value) {
  value <- gsub("&", "&amp;", as.character(value), fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  gsub(">", "&gt;", value, fixed = TRUE)
}

# Corresponds to generate_qc_report: Generate scan, schema, TIC and spike-in summary HTML.
generate_qc_report <- function(job) {
  mzml_file <- job[[1L]]; output_file <- job[[2L]]; spikeins <- job[[3L]]
  if (is.null(spikeins)) {
    cat("Using default spike-ins for ", as.character(mzml_file), "\n", sep = "")
    spikeins <- list(
      list("13C6-D-glucos", 187.0908, 0), list("trimethyl-13C3-caffeine", 198.0977, 0),
      list("15N-13C5-methionine", 156.0721, 0), list("13C5-L-glutamate", 153.0722, 0),
      list("15N2-uracil", 115.0286, 0), list("15N-L-tyrosine", 183.0782, 0)
    )
  } else if (is.character(spikeins) && length(spikeins) == 1L && endsWith(spikeins, "json")) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading spike-in JSON requires jsonlite.")
    spikeins <- jsonlite::fromJSON(spikeins, simplifyVector = FALSE)
  } else if (!is.list(spikeins)) {
    stop("Spike-ins must be a list of tuples or a JSON file.")
  }
  spectra <- .ms2_read_spectra(mzml_file)
  levels <- table(vapply(spectra, `[[`, 0L, "ms_level"))
  modes <- table(vapply(spectra, function(spectrum) {
    polarity <- .qc_or(spectrum$polarity, 0L)
    if (isTRUE(polarity == 1L) || identical(spectrum$positive_mode, TRUE)) "Positive" else "Negative"
  }, ""))
  description <- paste0(
    "The single file QC is a quick screening tool to check issues in mzML files.<br>",
    "This includes: the TIC, histograms of intensity and mz, and examination of spike-in standards.<br>",
    "Scan summary: ", paste(paste0(levels, " scans for ms_level=", names(levels)), collapse = ","), "<br>",
    "Mode summary: ", paste(paste0(modes, " scans for ion mode=", names(modes)), collapse = ","), "<br>"
  )
  tic <- .qc_calcTIC(spectra)
  trios <- .qc_extract_trio(spectra)
  solutions <- .qc_find_targets(trios, spikeins)
  rows <- character()
  for (name in names(solutions)) {
    target <- spikeins[[which(vapply(spikeins, `[[`, "", 1L) == name)[[1L]]]]
    found <- solutions[[name]]
    delta <- if (length(found)) vapply(found, `[[`, 0, 2L) - target[[2L]] else numeric()
    intensities <- if (length(found)) vapply(found, `[[`, 0, 1L) else numeric()
    best <- if (length(found)) found[[which.max(intensities)]] else c(NA, NA, NA)
    ppm <- if (length(delta)) sum(delta) / length(delta) / target[[2L]] * 1e6 else "NA"
    # The original version of Python writes the second item m/z of best tuple to max_intensity_time, retaining this behavior.
    rows <- c(rows, paste0(
      "<tr><td>", .qc_html(name), "</td><td>", ifelse(length(found), best[[1L]], "NA"),
      "</td><td>", ifelse(length(found), best[[2L]], "NA"), "</td><td>", ppm,
      "</td><td>", length(found) > 0L, "</td></tr>"
    ))
  }
  tic_text <- paste(paste(tic[[1L]], tic[[2L]], sep = ":"), collapse = ", ")
  html <- paste0(
    "<!DOCTYPE html><html lang='en'><head><meta charset='UTF-8'><title>Single File QC for ",
    .qc_html(mzml_file), "</title></head><body><h1>Single File QC for ", .qc_html(mzml_file),
    "</h1><p>", description, "</p><h2>Total Ion Chromatogram</h2><p>", tic_text,
    "</p><table><thead><tr><th>target</th><th>max_intensity</th><th>max_intensity_time</th>",
    "<th>ppm_error</th><th>detected</th></tr></thead><tbody>", paste(rows, collapse = ""),
    "</tbody></table></body></html>"
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  output_file
}

# Corresponds to get_dataframe_from_file: index_col and header retain Python0 base parameter semantics.
get_dataframe_from_file <- function(
    infile, header = 0L, index_col = 0L, sep = "\t", max_col = 21L) {
  data <- utils::read.table(
    infile, header = identical(as.integer(header), 0L), sep = sep,
    row.names = as.integer(index_col) + 1L, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(data) > max_col) {
    cat(sprintf("[asari.qc] Warning: too many columns (%d), truncating to %d.\n", ncol(data), max_col))
    data <- data[, seq_len(max_col), drop = FALSE]
  }
  data
}
