# Corresponds to Python asari/tools/plot_scan_seq.py: MS1/MS2 scan sequence summary plot.

# Corresponds to summarize_ms_file.
summarize_ms_file <- function(
    infile, nspec_plot = 50L, offset_n = 100L, width = 10, height = 6,
    title = "", outfile = "ms_summary.pdf") {
  spectra <- extract_all_spectra_form_file(infile, min_intensity = 100, MS2_peak_limit = 999999L)
  plot_msn(
    spectra[[1L]], spectra[[2L]], nspec_plot, offset_n,
    width, height, title, outfile
  )
  invisible(NULL)
}

# Corresponds to plot_msn internal get_median_height.
.plot_msn_get_median_height <- function(peaks) {
  # Empty scans have no definable median intensity and return NA for skipping by the drawing stage.
  if (length(peaks) == 0L) return(NA_real_)
  log10(stats::median(vapply(peaks, `[[`, 0, 2L)))
}

# Corresponds to plot_msn internal get_max_height.
.plot_msn_get_max_height <- function(peaks) {
  if (length(peaks)) max(vapply(peaks, `[[`, 0, 2L)) else 0
}

# Corresponds to plot_msn: number of scanned peaks/median intensity in the upper figure, MS1 base peak XIC in the lower figure.
plot_msn <- function(
    ms1_spectra, ms2_spectra, nspec_plot = 50L, offset_n = 100L,
    width = 10, height = 6, title = "", outfile = NULL) {
  number1 <- length(ms1_spectra); number2 <- length(ms2_spectra)
  if (number1 == 0L) stop("No MS1 spectrum is available for the scan summary.", call. = FALSE)
  all_spectra <- c(ms1_spectra, ms2_spectra)
  max_rtime <- max(vapply(all_spectra, `[[`, 0, "rtime"))
  indices <- if (number1 + number2 < nspec_plot + offset_n) 0:(nspec_plot - 1L) else offset_n:(offset_n + nspec_plot - 1L)
  ids <- paste0("sp", indices + 1L)
  decorate <- function(spectra, level) lapply(Filter(function(item) item$id %in% ids, spectra), function(item) {
    item$level <- level; item$num_peaks <- length(item$peaks)
    item$median_intensity <- .plot_msn_get_median_height(item$peaks); item
  })
  plot1 <- decorate(ms1_spectra, 1L); plot2 <- decorate(ms2_spectra, 2L)
  if (!is.null(outfile)) grDevices::pdf(outfile, width = width, height = height)
  if (!is.null(outfile)) on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(mfrow = c(2, 1), oma = c(0, 0, 2, 0)); on.exit(graphics::par(old), add = TRUE)
  all_plot <- c(plot1, plot2)
  # When the offset exceeds the short file, it returns to the front scan, ensuring that small or pure MS1 files can also be output.
  if (length(all_plot) == 0L) {
    ids <- vapply(utils::head(all_spectra, nspec_plot), `[[`, "", "id")
    plot1 <- decorate(ms1_spectra, 1L); plot2 <- decorate(ms2_spectra, 2L)
    all_plot <- c(plot1, plot2)
  }
  finite_plot <- Filter(function(item) is.finite(item$median_intensity), all_plot)
  if (length(finite_plot) == 0L) {
    stop("Selected scans do not contain peaks above the plotting threshold.", call. = FALSE)
  }
  median_values <- vapply(finite_plot, `[[`, 0, "median_intensity")
  y_lower <- min(2, min(median_values))
  y_max <- max(median_values) + 0.5
  graphics::plot(NA, xlim = range(vapply(finite_plot, `[[`, 0, "rtime")), ylim = c(y_lower, y_max), xlab = "Retention time (sec)", ylab = "Median log10 intensity", main = "Example scans")
  for (item in plot1) graphics::segments(item$rtime, 2, item$rtime, item$median_intensity, col = "blue")
  for (item in plot2) graphics::segments(item$rtime, 2, item$rtime, item$median_intensity, col = "red", lty = 2)
  graphics::plot(
    vapply(ms1_spectra, `[[`, 0, "rtime"),
    vapply(ms1_spectra, function(item) .plot_msn_get_max_height(item$peaks), 0),
    type = "l", col = "blue", xlab = "Retention time (sec)", ylab = "Base peak intensity", main = "Ion chromatogram"
  )
  graphics::legend("topleft", legend = c(
    paste("Number of MS1 scans:", number1), paste("Number of MS2 scans:", number2),
    sprintf("Max retention time: %.2f sec", max_rtime)
  ), bty = "n")
  graphics::mtext(title, outer = TRUE, side = 3, adj = 0.05)
  invisible(NULL)
}
