# Corresponds to Python asari/tools/plot.py: mass-spectral regions, mass track, mirror spectrum and QC plot.

# Corresponds to get_plot_region_from_file: both scan number and m/z use strict open intervals.
get_plot_region_from_file <- function(
    infile, min_scan_number, max_scan_number, min_mz, max_mz, ms_level = 1L) {
  spectra <- .ms2_read_spectra(infile)
  result <- list()
  for (ii in seq_along(spectra)) {
    scan0 <- ii - 1L
    spectrum <- spectra[[ii]]
    if (min_scan_number < scan0 && scan0 < max_scan_number && spectrum$ms_level == ms_level) {
      peaks <- .ms2_peak_matrix(spectrum$peaks)
      keep <- peaks[, 1L] > min_mz & peaks[, 1L] < max_mz
      for (row in which(keep)) result[[length(result) + 1L]] <- c(scan0, peaks[row, 1L], as.integer(peaks[row, 2L]))
    }
  }
  result
}

# Convert the datapoint list into a three-column matrix and calculate the relative log2 intensity.
.plot_region_data <- function(datapoints) {
  data <- do.call(rbind, lapply(datapoints, as.numeric))
  intensity <- log2(data[, 3L] + 1)
  normalized <- intensity / max(intensity)
  list(data = data, normalized = normalized)
}

# Corresponds to plot_scatter_map_region.
plot_scatter_map_region <- function(datapoints, figsize = c(8, 10), cmap = NULL,
                                    colorbar_orientation = "horizontal") {
  invisible(figsize); invisible(cmap); invisible(colorbar_orientation)
  prepared <- .plot_region_data(datapoints)
  graphics::plot(
    prepared$data[, 1L], prepared$data[, 2L], pch = 20,
    col = grDevices::hcl.colors(100, "Blue-Red 3")[pmax(1L, ceiling(prepared$normalized * 100))],
    xlab = "scan number", ylab = "m/z"
  )
  invisible(NULL)
}

# Corresponds to double_scatter_map_region.
double_scatter_map_region <- function(datapoints, figsize = c(8, 10), cmap = NULL,
                                      colorbar_orientation = "horizontal") {
  invisible(figsize); invisible(cmap); invisible(colorbar_orientation)
  prepared <- .plot_region_data(datapoints)
  old <- graphics::par(mfrow = c(2, 1)); on.exit(graphics::par(old), add = TRUE)
  colors <- grDevices::hcl.colors(100, "Blue-Red 3")[pmax(1L, ceiling(prepared$normalized * 100))]
  graphics::plot(prepared$data[, 1L], prepared$data[, 2L], pch = 20, col = colors, xlab = "", ylab = "m/z")
  graphics::plot(prepared$data[, 1L], prepared$data[, 3L], pch = 20, col = colors, xlab = "scan number", ylab = "intensity")
  invisible(NULL)
}

# Corresponds to with_line_scatter_map_region.
with_line_scatter_map_region <- function(datapoints, figsize = c(8, 10), cmap = NULL) {
  invisible(figsize); invisible(cmap)
  prepared <- .plot_region_data(datapoints)
  old <- graphics::par(mfrow = c(2, 1)); on.exit(graphics::par(old), add = TRUE)
  graphics::plot(prepared$data[, 1L], prepared$data[, 2L], pch = 20, xlab = "", ylab = "m/z")
  graphics::plot(prepared$data[, 1L], prepared$data[, 3L], type = "o", pch = 8, xlab = "scan number", ylab = "intensity")
  invisible(NULL)
}

# Corresponds to plot_masstrack: start/end follows Python half-open 0-based slicing.
plot_masstrack <- function(track, color = "magenta", start = 100L, end = 400L,
                           yticks = c(0, 5e7, 1e8)) {
  x <- if (end > start) start:(end - 1L) else integer()
  y <- track$intensity[x + 1L]
  graphics::plot(x, y, col = color, type = "o", pch = 1, lty = 2, yaxt = "n", xlab = "scan", ylab = "intensity")
  graphics::axis(2, at = yticks)
  invisible(NULL)
}

# Map matplotlib color names to R-recognizable colors.
.plot_color <- function(value) {
  switch(value, "tab:red" = "red3", "R0" = "red", value)
}

# Corresponds to mirror_plot: normalizes, flips up and down, and joins within-tolerance peaks.
mirror_plot <- function(
    peaks, peaks2, figsize = c(8, 4), label1 = "GCMS features",
    label2 = "Spectrum in lib", normalize = TRUE, match_tol = NULL,
    colors = c("blue", "tab:red", "black"), title = "GCMS Mirror Plot",
    outfile = "this_mirror_plot.pdf") {
  first <- if (is.matrix(peaks)) peaks else do.call(rbind, peaks)
  second <- if (is.matrix(peaks2)) peaks2 else do.call(rbind, peaks2)
  i1 <- first[, 2L]; i2 <- second[, 2L]
  if (isTRUE(normalize)) {
    if (max(i1) > 0) i1 <- i1 / max(i1)
    if (max(i2) > 0) i2 <- i2 / max(i2)
  }
  opened <- !is.null(outfile) && nzchar(outfile)
  if (opened) grDevices::pdf(outfile, width = figsize[[1L]], height = figsize[[2L]])
  if (opened) on.exit(grDevices::dev.off(), add = TRUE)
  xlim <- range(c(first[, 1L], second[, 1L]))
  graphics::plot(NA, xlim = xlim, ylim = c(-1.05, 1.05), xlab = "m/z", ylab = "Normalized Intensity", main = title)
  graphics::segments(first[, 1L], 0, first[, 1L], i1, col = .plot_color(colors[[1L]]))
  graphics::segments(second[, 1L], 0, second[, 1L], -i2, col = .plot_color(colors[[2L]]))
  if (!is.null(match_tol)) {
    for (ii in seq_len(nrow(first))) for (jj in which(abs(second[, 1L] - first[ii, 1L]) <= match_tol)) {
      graphics::segments(first[ii, 1L], i1[[ii]], second[jj, 1L], -i2[[jj]], col = .plot_color(colors[[3L]]), lwd = 0.3)
    }
  }
  graphics::abline(h = 0)
  graphics::legend("topleft", legend = c(label1, label2), col = vapply(colors[1:2], .plot_color, ""), lty = 1, bty = "n")
  invisible(NULL)
}

# Corresponds to plot_mSelectivity.
plot_mSelectivity <- function(mzList, selectivities, figsize = c(10, 3), save_pdf = TRUE,
                              outfile = "Figure_mSelectivity") {
  if (isTRUE(save_pdf)) grDevices::pdf(paste0(outfile, ".pdf"), width = figsize[[1L]], height = figsize[[2L]])
  if (isTRUE(save_pdf)) on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(mzList, selectivities, type = "o", xlab = "m/z", ylab = "Selectivity", main = "m/z selectivity (zoom in)")
  graphics::segments(mzList, 0, mzList, 1, col = "red", lwd = 0.2)
  invisible(NULL)
}

# Corresponds to plot_cSelectivity.
plot_cSelectivity <- function(mass_track, list_ranges = list(c(0, 800), c(100, 150)),
                              save_pdf = TRUE, outfile = "Figure_cSelectivity") {
  if (isTRUE(save_pdf)) grDevices::pdf(paste0(outfile, ".pdf"), width = 10, height = 4)
  if (isTRUE(save_pdf)) on.exit(grDevices::dev.off(), add = TRUE)
  full <- list_ranges[[1L]]
  x <- if (full[[2L]] > full[[1L]]) full[[1L]]:(full[[2L]] - 1L) else integer()
  graphics::plot(x, mass_track$intensity[x + 1L], pch = 20, col = "black", xlab = "scan", ylab = "intensity")
  if (length(list_ranges) > 1L) for (range in list_ranges[-1L]) {
    x2 <- range[[1L]]:(range[[2L]] - 1L)
    graphics::polygon(c(x2, rev(x2)), c(rep(0, length(x2)), rev(mass_track$intensity[x2 + 1L])), col = grDevices::adjustcolor("red", 0.2), border = NA)
  }
  invisible(NULL)
}

# Corresponds to pca_ftable: Calculate the two principal components after normalizing the sample transpose matrix.
pca_ftable <- function(tsv_path) {
  table <- utils::read.delim(tsv_path, check.names = FALSE)
  samples <- t(as.matrix(table[, 12:ncol(table), drop = FALSE]))
  pca <- stats::prcomp(samples, center = TRUE, scale. = TRUE, rank. = 2L)
  scores <- pca$x[, 1:2, drop = FALSE]
  graphics::plot(scores[, 1L], scores[, 2L], type = "n", xlab = "PC1", ylab = "PC2")
  graphics::text(scores[, 1L], scores[, 2L], labels = rownames(scores))
  invisible(scores)
}

# Corresponds to plot_correlations: Compare the sample aggregation intensity with the comprehensive peak_area.
plot_correlations <- function(feature_table_path) {
  table <- utils::read.delim(feature_table_path, check.names = FALSE)
  samples <- log2(as.matrix(table[, 12:ncol(table), drop = FALSE]) + 1)
  summaries <- list(
    max = apply(samples, 1L, max), mean = rowMeans(samples),
    median = apply(samples, 1L, stats::median), min = apply(samples, 1L, min)
  )
  peak_area <- log2(table$peak_area + 1)
  old <- graphics::par(mfrow = c(2, 2)); on.exit(graphics::par(old), add = TRUE)
  for (name in names(summaries)) {
    cat(name, "\n")
    x <- summaries[[name]]
    r2 <- round(stats::cor(x, peak_area)^2, 2)
    graphics::plot(x, peak_area, xlab = name, ylab = "peak_area", main = paste0(name, " vs peak_area (R2 = ", r2, ")"))
    graphics::abline(stats::lm(peak_area ~ x), col = "blue")
  }
  invisible(NULL)
}

# Read object fields, compatible with environment/list.
.plot_field <- function(object, name) if (is.environment(object)) object[[name]] else object[[name]]

# Corresponds to plot_peaks_mastrace: plots the raw mass track and fitted peak model.
plot_peaks_masstrace <- function(sample, mzstr, outfile = "masstrace_plot.pdf") {
  grDevices::pdf(outfile); on.exit(grDevices::dev.off(), add = TRUE)
  traces <- .plot_field(sample, "dict_masstraces")[[mzstr]]
  if (length(traces) == 0L) stop("No mass traces found for ", mzstr)
  first_trace <- traces[[1L]]
  graphics::plot(
    .plot_field(first_trace, "list_retention_time"),
    .plot_field(first_trace, "list_intensity"),
    pch = 20, cex = 0.3, xlab = "retention time", ylab = "intensity"
  )
  if (length(traces) > 1L) for (trace in traces[-1L]) {
    graphics::points(
      .plot_field(trace, "list_retention_time"),
      .plot_field(trace, "list_intensity"), pch = 20, cex = 0.3
    )
  }
  creator <- .plot_field(sample, "create_peak_dict")
  peaks <- creator()[[mzstr]]
  for (peak in peaks) {
    extender <- .plot_field(peak, "extend_model_range")
    extender()
    graphics::lines(.plot_field(peak, "rt_extended"), .plot_field(peak, "y_fitted_extended"), col = grDevices::adjustcolor("red", 0.5), lwd = 0.6)
  }
  graphics::title(paste("mass trace", round(.plot_field(tail(traces, 1L)[[1L]], "mz"), 6)))
  invisible(NULL)
}

# Corresponds to the current Python empty placeholder plot_peaks.
plot_peaks <- function() invisible(NULL)

# Corresponds to plot_sample_rt_calibration.
plot_sample_rt_calibration <- function(sample, outfile = "rt_calibration.pdf") {
  path <- paste0(.plot_field(sample, "name"), outfile)
  grDevices::pdf(path); on.exit(grDevices::dev.off(), add = TRUE)
  data <- .plot_field(sample, "__rt_calibration__data__")
  graphics::plot(data[[1L]], data[[2L]], type = "o", pch = 1, cex = 0.3, main = "rt_calibration")
  invisible(NULL)
}
