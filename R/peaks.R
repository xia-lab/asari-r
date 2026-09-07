# Corresponds to Python asari/peaks.py: chromatographic peak detection, noise estimation and peak quality evaluation.
#
# Compatibility conventions:
# 1. The externally exposed scan/apex/left_base/right_base all maintain Python’s 0-based numbering.
# 2. Mass track uses a named list, and the fields are consistent with Python dict.
# 3. This file uses base R and stats to avoid introducing large dependencies for a small number of SciPy functions.

# Replicate the behavior of Python slicing [start:stop] in R.
# Both start and stop use Python's 0-based bounds; the position of stop is not included in the result.
.python_slice <- function(values, start = NULL, stop = NULL) {
  n <- length(values)

  # Python interprets negative bounds as counting backwards from the end of the vector.
  normalize_bound <- function(bound, default) {
    if (is.null(bound)) return(default)
    bound <- as.integer(bound)
    if (bound < 0L) bound <- max(n + bound, 0L)
    min(max(bound, 0L), n)
  }

  start <- normalize_bound(start, 0L)
  stop <- normalize_bound(stop, n)
  if (stop <= start) return(values[0])

  # R subscripts start at 1, so convert Python's [start, stop) to R intervals.
  values[seq.int(start + 1L, stop)]
}

# Reproduce Python int()'s truncation toward zero while avoiding R 32-bit integer overflow NA.
.python_int <- function(value) {
  truncated <- trunc(value)
  if (length(truncated) == 1L && is.finite(truncated) &&
      truncated >= -.Machine$integer.max &&
      truncated <= .Machine$integer.max) {
    return(as.integer(truncated))
  }

  # When the range of R integer is exceeded, use non-decimal double to save the same integer value.
  truncated
}

# Reproduce linear detrending of scipy.signal.detrend(..., type = "linear").
.detrend_linear <- function(values) {
  values <- as.numeric(values)
  if (length(values) < 2L) return(values - mean(values))

  # SciPy fits a straight line to the sample number, and then returns the original signal minus the fitted value.
  x <- seq_along(values) - 1L
  fit <- stats::lm.fit(cbind(1, x), values)
  values - as.vector(cbind(1, x) %*% fit$coefficients)
}

# Find local maxima and plateau peak centers, the behavior corresponds to the first stage of SciPy find_peaks.
.local_maxima_1d <- function(values) {
  values <- as.numeric(values)
  n <- length(values)
  if (n < 3L) return(integer())

  peaks <- integer()
  ii <- 2L
  while (ii <= n - 1L) {
    # The left side of the peak must strictly rise; equal values are handled uniformly by the subsequent platform logic.
    if (values[[ii - 1L]] < values[[ii]]) {
      right_edge <- ii
      while (right_edge < n && values[[right_edge + 1L]] == values[[ii]]) {
        right_edge <- right_edge + 1L
      }

      # The right side of the peak must also go down strictly; plateau peaks take the midpoint and round down.
      if (right_edge < n && values[[right_edge]] > values[[right_edge + 1L]]) {
        peaks <- c(peaks, as.integer(floor((ii + right_edge) / 2)))
      }
      ii <- right_edge + 1L
    } else {
      ii <- ii + 1L
    }
  }

  # The private peak finder internally uses R's 1-based subscripting.
  peaks
}

# High priority peaks are retained according to SciPy's distance rule.
.select_by_peak_distance <- function(peaks, heights, distance) {
  if (length(peaks) < 2L) return(seq_along(peaks))
  distance <- as.integer(ceiling(distance))
  keep <- rep(TRUE, length(peaks))

  # SciPy prioritizes peak heights; the same peaks are processed from right to left after being sorted stably by position.
  priority_order <- order(heights, seq_along(heights))
  for (position in rev(priority_order)) {
    if (!keep[[position]]) next

    # Delete left and right neighboring peaks that are less than distance from the current high-priority peak.
    left <- position - 1L
    while (left >= 1L && peaks[[position]] - peaks[[left]] < distance) {
      keep[[left]] <- FALSE
      left <- left - 1L
    }
    right <- position + 1L
    while (right <= length(peaks) && peaks[[right]] - peaks[[position]] < distance) {
      keep[[right]] <- FALSE
      right <- right + 1L
    }
  }

  which(keep)
}

# Calculate prominence, left and right base, and half-maximum width for a single peak.
.peak_prominence_and_width <- function(values, peak, wlen = NULL) {
  n <- length(values)
  peak_height <- values[[peak]]

  # SciPy will round wlen up to an odd number, and center the window around the peak.
  if (is.null(wlen)) {
    left_limit <- 1L
    right_limit <- n
  } else {
    rounded_wlen <- as.integer(ceiling(wlen))
    if (rounded_wlen %% 2L == 0L) rounded_wlen <- rounded_wlen + 1L
    half_window <- rounded_wlen %/% 2L
    left_limit <- max(1L, peak - half_window)
    right_limit <- min(n, peak + half_window)
  }

  # Scan left from the peak; stop at a higher point and record the lowest point as left base.
  left_base <- peak
  left_min <- peak_height
  cursor <- peak
  while (cursor >= left_limit && values[[cursor]] <= peak_height) {
    if (values[[cursor]] < left_min) {
      left_min <- values[[cursor]]
      left_base <- cursor
    }
    cursor <- cursor - 1L
  }

  # Perform the same scan from the peak to the right to get the right base.
  right_base <- peak
  right_min <- peak_height
  cursor <- peak
  while (cursor <= right_limit && values[[cursor]] <= peak_height) {
    if (values[[cursor]] < right_min) {
      right_min <- values[[cursor]]
      right_base <- cursor
    }
    cursor <- cursor + 1L
  }

  prominence <- peak_height - max(left_min, right_min)
  evaluation_height <- peak_height - 0.5 * prominence

  # Finds the half-height intersection between the left base and the peak and linearly interpolates.
  left_ip <- as.numeric(peak)
  cursor <- peak
  while (cursor > left_base && values[[cursor]] > evaluation_height) {
    cursor <- cursor - 1L
  }
  if (values[[cursor]] < evaluation_height && cursor < peak) {
    left_ip <- cursor + (evaluation_height - values[[cursor]]) /
      (values[[cursor + 1L]] - values[[cursor]])
  } else {
    left_ip <- as.numeric(cursor)
  }

  # Finds the half-height intersection between the peak and right base and linearly interpolates.
  right_ip <- as.numeric(peak)
  cursor <- peak
  while (cursor < right_base && values[[cursor]] > evaluation_height) {
    cursor <- cursor + 1L
  }
  if (values[[cursor]] < evaluation_height && cursor > peak) {
    right_ip <- cursor - (evaluation_height - values[[cursor]]) /
      (values[[cursor - 1L]] - values[[cursor]])
  } else {
    right_ip <- as.numeric(cursor)
  }

  list(
    prominence = prominence,
    left_base = left_base,
    right_base = right_base,
    width = right_ip - left_ip,
    left_ip = left_ip,
    right_ip = right_ip
  )
}

# Implement the subset of scipy.signal.find_peaks parameters actually used by peaks.py.
.find_peaks_scipy <- function(values,
                              height = NULL,
                              distance = NULL,
                              prominence = NULL,
                              width = NULL,
                              wlen = NULL) {
  values <- as.numeric(values)

  # The first step is to find all local maxima and plateau peak centers.
  peaks <- .local_maxima_1d(values)
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(), prominences = numeric(),
      left_bases = integer(), right_bases = integer(), widths = numeric()
    )))
  }

  # The second step applies the minimum peak height; here the 1-based peak positions within R are still retained.
  if (!is.null(height)) {
    peaks <- peaks[values[peaks] >= height]
  }
  if (length(peaks) == 0L) {
    return(list(peaks = integer(), properties = list(
      peak_heights = numeric(), prominences = numeric(),
      left_bases = integer(), right_bases = integer(), widths = numeric()
    )))
  }

  # The third step is to perform minimum peak spacing filtering based on peak height priority.
  if (!is.null(distance)) {
    selected <- .select_by_peak_distance(peaks, values[peaks], distance)
    peaks <- peaks[selected]
  }

  # The fourth step is to calculate prominence, bases and FWHM width.
  metrics <- lapply(peaks, function(peak) {
    .peak_prominence_and_width(values, peak, wlen)
  })
  prominences <- vapply(metrics, `[[`, numeric(1), "prominence")
  widths <- vapply(metrics, `[[`, numeric(1), "width")

  # SciPy filters prominence first, then width; it filters all properties simultaneously each time.
  keep <- rep(TRUE, length(peaks))
  if (!is.null(prominence)) keep <- keep & prominences >= prominence
  if (!is.null(width)) keep <- keep & widths >= width
  peaks <- peaks[keep]
  metrics <- metrics[keep]
  prominences <- prominences[keep]
  widths <- widths[keep]

  # External conversion back to 0-based peak and base numbers used by Python.
  properties <- list(
    peak_heights = values[peaks],
    prominences = prominences,
    left_bases = vapply(metrics, `[[`, integer(1), "left_base") - 1L,
    right_bases = vapply(metrics, `[[`, integer(1), "right_base") - 1L,
    widths = widths,
    left_ips = vapply(metrics, `[[`, numeric(1), "left_ip") - 1,
    right_ips = vapply(metrics, `[[`, numeric(1), "right_ip") - 1
  )

  list(peaks = as.integer(peaks - 1L), properties = properties)
}

# Fits a Gaussian using a nonlinear least squares objective; on failure the same Python fallback value is taken by the caller.
.fit_gaussian <- function(x, y, start) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  start <- as.numeric(start)

  # The goal of curve_fit is to minimize the sum of squared residuals between Gaussian predictions and observations.
  objective <- function(parameters) {
    predicted <- gaussian_function__(
      x, parameters[[1L]], parameters[[2L]], parameters[[3L]]
    )
    residual_sum <- sum((predicted - y)^2)
    if (is.finite(residual_sum)) residual_sum else .Machine$double.xmax
  }

  # Take the original fast BFGS path first; most peaks can get stable solutions here.
  parameter_scale <- pmax(abs(start), c(1, 1, 1))
  fit <- stats::optim(
    start,
    objective,
    method = "BFGS",
    control = list(
      maxit = 10000L,
      reltol = 1e-12,
      parscale = parameter_scale
    )
  )
  fit_error <- objective(fit$par)
  centered_error <- sum((y - mean(y))^2)
  fit_is_valid <- fit$convergence == 0L && all(is.finite(fit$par)) &&
    fit_error <= objective(start) && fit_error <= centered_error &&
    abs(fit$par[[3L]]) <= 10 * max(1, diff(range(x)))

  if (fit_is_valid) {
    # Perform a small number of Gauss-Newton iterations from the BFGS solution to get close to the machine accuracy results of curve_fit.
    refined <- tryCatch(
      suppressWarnings(stats::nls(
        y ~ a * exp(-((x - mu)^2) / (2 * sigma^2)),
        start = list(a = fit$par[[1L]], mu = fit$par[[2L]], sigma = fit$par[[3L]]),
        control = stats::nls.control(maxiter = 20L, tol = 1e-10, warnOnly = TRUE)
      )),
      error = function(e) NULL
    )
    if (!is.null(refined) && all(is.finite(stats::coef(refined)))) {
      return(stats::coef(refined))
    }
    return(stats::setNames(fit$par, c("a", "mu", "sigma")))
  }

  # Robust recovery is enabled from the same initial values in SciPy only if the fast solution does not improve the initial values.
  direct_nls <- tryCatch(
    suppressWarnings(stats::nls(
      y ~ a * exp(-((x - mu)^2) / (2 * sigma^2)),
      start = list(a = start[[1L]], mu = start[[2L]], sigma = start[[3L]]),
      control = stats::nls.control(
        maxiter = 1000L, tol = 1e-12, minFactor = 1e-8, warnOnly = TRUE
      )
    )),
    error = function(e) NULL
  )
  if (!is.null(direct_nls) && all(is.finite(stats::coef(direct_nls)))) {
    direct_parameters <- unname(stats::coef(direct_nls))
    if (objective(direct_parameters) <= objective(start)) {
      return(stats::setNames(direct_parameters, c("a", "mu", "sigma")))
    }
  }

  stop("Gaussian fitting failed to converge")
}

# Detect multiple mass tracks in batches and flatten the results of each track into a peak list.
batch_deep_detect_elution_peaks <- function(list_mass_tracks,
                                            number_of_scans,
                                            parameters) {
  # First generate a job list consistent with the Python multiprocessing input.
  jobs <- iter_peak_detection_parameters(
    list_mass_tracks, number_of_scans, parameters
  )

  # peak detection tasks are independent of each other. After being parallelized by multicores, the input order is still maintained by bulk_process.
  # Low-level direct calls that lack this parameter continue to use a single core, avoiding changes to existing testing and calling behavior.
  workers <- parameters$multicores
  if (is.null(workers) || length(workers) != 1L || is.na(workers)) workers <- 1L
  workers <- max(1L, min(as.integer(workers), length(jobs)))
  results <- bulk_process(
    stats_detect_elution_peaks,
    jobs,
    jobs_per_worker = workers
  )
  if (length(results) == 0L) return(list())
  unlist(results, recursive = FALSE)
}

# Generate peak detection task for mass tracks that exceed a minimum intensity threshold.
iter_peak_detection_parameters <- function(list_mass_tracks,
                                           number_of_scans,
                                           parameters) {
  # Python np.any(...) retains the track as long as any intensity is strictly greater than the threshold.
  selected_tracks <- Filter(function(track) {
    any(track$intensity > parameters$min_intensity_threshold)
  }, list_mass_tracks)

  # Each job maintains a (mass_track, number_of_scans, parameters) ternary structure.
  lapply(selected_tracks, function(track) {
    list(track, number_of_scans, parameters)
  })
}

# Execute asari's statistical guide peak detection main process for a mass track.
stats_detect_elution_peaks <- function(job) {
  # Unpack the Python job triplet and read all detection parameters.
  mass_track <- job[[1L]]
  number_of_scans <- as.integer(job[[2L]])
  parameters <- job[[3L]]
  min_peak_height <- parameters$min_peak_height
  min_peak_ratio <- parameters$min_peak_ratio
  min_fwhm <- round(0.5 * parameters$min_timepoints)
  min_intensity_threshold <- parameters$min_intensity_threshold
  wlen <- parameters$wlen
  snr <- parameters$signal_noise_ratio
  peakshape <- parameters$gaussian_shape
  min_prominence_ratio <- 0.02
  min_prominence_threshold <- parameters$min_prominence_threshold

  # First clean, scale and estimate the baseline and noise of the entire mass track.
  audited <- audit_mass_track(
    mass_track$intensity,
    min_fwhm,
    min_intensity_threshold,
    min_peak_height,
    min_peak_ratio
  )
  baseline <- audited[[1L]]
  noise_level <- audited[[2L]]
  scaling_factor <- audited[[3L]]
  min_peak_height <- audited[[4L]]
  list_intensity <- audited[[5L]]
  list_scans <- seq_len(number_of_scans) - 1L
  list_json_peaks <- list()
  list_peaks <- list()

  # Establish the ROI using scan points above the noise level; scan numbers remain 0-based.
  selected_scans <- list_scans[list_intensity > noise_level]

  # Python `array.any()` returns FALSE for arrays containing only zeros, and this behavior is explicitly reproduced here.
  if (length(selected_scans) > 0L && any(selected_scans != 0L)) {
    rois <- list()
    current_roi <- selected_scans[[1L]]

    # Adjacent signals allow for gaps in up to two scans and fill the gaps into the ROI.
    if (length(selected_scans) > 1L) {
      for (scan in selected_scans[-1L]) {
        if (scan - tail(current_roi, 1L) < 3L) {
          current_roi <- c(
            current_roi,
            seq.int(tail(current_roi, 1L) + 1L, scan)
          )
        } else {
          rois[[length(rois) + 1L]] <- current_roi
          current_roi <- scan
        }
      }
    }
    rois[[length(rois) + 1L]] <- current_roi

    # Delete ROIs that are too short to avoid peak width and Gaussian fitting without enough points.
    rois <- Filter(function(roi) length(roi) >= min_fwhm + 2L, rois)

    # Each ROI is detected separately and the peak records are summarized in JSON style.
    for (roi in rois) {
      if (length(roi) < 3 * min_fwhm) {
        roi <- extend_ROI(roi, number_of_scans)
      }
      intensity_roi <- list_intensity[roi + 1L]
      detected <- detect_evaluate_peaks_on_roi(
        intensity_roi,
        roi,
        min_peak_height,
        min_peak_ratio,
        min_fwhm,
        min_prominence_threshold,
        wlen,
        snr,
        peakshape,
        min_prominence_ratio,
        noise_level
      )
      list_json_peaks <- c(list_json_peaks, detected)
    }

    # The cSelectivity of the peaks is calculated and noisy scan points that do not belong to any peak are collected.
    selectivity_stats <- `__peaks_cSelectivity_stats_`(
      list_intensity, list_json_peaks
    )
    list_cSelectivity <- selectivity_stats[[1L]]
    peak_datapoints <- selectivity_stats[[2L]]
    noise_data_points <- setdiff(list_scans, peak_datapoints)

    # Write selectivity, parent mass track information, raw magnitude height, and local SNR back to each peak.
    for (ii in seq_along(list_json_peaks)) {
      peak <- list_json_peaks[[ii]]
      peak$cSelectivity <- list_cSelectivity[[ii]]
      peak$parent_masstrack_id <- mass_track$id_number
      peak$mz <- mass_track$mz
      peak$height <- scaling_factor * peak$height + baseline
      capped_height <- min(peak$height, 1e100)
      local_noise <- compute_noise_by_flanks(
        peak,
        mass_track$intensity,
        noise_data_points,
        min_intensity_threshold,
        noise_level
      )
      peak$snr <- .python_int(capped_height / local_noise)

      # Only peaks reaching the specified SNR are entered into the final result; heights are truncated as Python int().
      if (peak$snr >= snr) {
        peak$height <- .python_int(peak$height)
        list_peaks[[length(list_peaks) + 1L]] <- peak
      }
    }
  }

  list_peaks
}

# Local noise is estimated using non-peak data points on either side of the peak.
compute_noise_by_flanks <- function(peak,
                                    list_intensity,
                                    noise_data_points,
                                    min_intensity_threshold,
                                    old_noise_level) {
  padding <- 30L
  number_of_points <- 100L

  # Sort after adding apex and take the first matching position using Python list.index().
  points <- sort(c(noise_data_points, peak$apex))
  center <- match(peak$apex, points) - 1L

  # Strict reproduction of two Python slices, including semantics for negative starting points counting from the end of the list.
  left_points <- .python_slice(
    points,
    center - padding - number_of_points,
    center - padding
  )
  right_points <- .python_slice(
    points,
    center + padding,
    center + number_of_points + padding
  )
  use_index <- c(left_points, right_points)

  # Python adds 1 to the intensity and takes the mean, using the lowest intensity threshold as the lower bound.
  if (length(use_index) > 0L) {
    return(max(
      min_intensity_threshold,
      mean(list_intensity[use_index + 1L] + 1)
    ))
  }
  old_noise_level
}

# Cleans the mass track and returns the baseline, noise, scaling factor, dynamic peak height threshold, and cleaned intensity.
audit_mass_track <- function(list_intensity,
                             min_fwhm,
                             min_intensity_threshold,
                             min_peak_height,
                             min_peak_ratio) {
  # SciPy uniform_filter1d will retain np.int64; R uses numeric to save traces to avoid 32-bit integer overflow.
  # So as long as the full intensity remains a finite integer value, Python's truncation semantics for integer output continue to apply.
  working_is_integer <- is.integer(list_intensity) || (
    is.numeric(list_intensity) && all(is.finite(list_intensity)) &&
      all(list_intensity == trunc(list_intensity))
  )
  list_intensity <- as.numeric(list_intensity)
  scaling_factor <- 1
  low <- min_intensity_threshold
  high <- 1e8
  baseline <- low
  noise_level <- low
  max_intensity <- max(list_intensity)

  # Tracks exceeding 1e8 are scaled first to prevent prominence calculations from being destroyed by extremely large numbers.
  if (max_intensity > high) {
    scaling_factor <- max_intensity / high
    list_intensity <- list_intensity / scaling_factor
    working_is_integer <- FALSE
  }

  # Baseline and noise are re-estimated only if the median exceeds a minimum threshold.
  median_intensity <- stats::median(list_intensity)
  if (median_intensity > low) {
    above_low <- list_intensity[list_intensity > low]

    # The high-baseline track first performs the same linear detrending as scipy.signal.detrend.
    if (length(above_low) > length(list_intensity) * 0.5 &&
        median_intensity > 10 * min_peak_height) {
      list_intensity <- .detrend_linear(list_intensity)
      working_is_integer <- FALSE
    }

    # Estimate the background distribution by taking points below the "lower quartile + LOW".
    bottom_cutoff <- low + as.numeric(stats::quantile(
      list_intensity, 0.25, names = FALSE, type = 7
    ))
    bottom_values <- list_intensity[list_intensity < bottom_cutoff]
    baseline <- mean(bottom_values)

    # NumPy std defaults to ddof=0, and R's sd(), which defaults to ddof=1, cannot be used directly.
    noise_level <- sqrt(mean((bottom_values - baseline)^2))
  }

  # Neither baseline nor noise can fall below a preset minimum intensity threshold.
  baseline <- max(baseline, low)
  noise_level <- max(noise_level, low)

  # Tracks that are extremely tall, scaled, or relatively noisy are moving averaged according to Python rules.
  if (scaling_factor > 1 || max_intensity > 1000 * min_peak_height) {
    list_intensity <- smooth_moving_average(
      list_intensity, size = min_fwhm * 3
    )
    if (working_is_integer) list_intensity <- trunc(list_intensity)
  } else if (100 * noise_level > max_intensity ||
             max_intensity < 10 * min_peak_height) {
    list_intensity <- smooth_moving_average(
      list_intensity, size = min_fwhm + 2
    )
    if (working_is_integer) list_intensity <- trunc(list_intensity)
  }

  # Finally, the baseline is subtracted, and the lowest peak height is dynamically adjusted according to the maximum intensity of the full track.
  list_intensity <- list_intensity - baseline
  min_peak_height <- max(
    min_peak_height,
    max_intensity * min_peak_ratio
  ) - baseline

  list(
    baseline,
    noise_level,
    scaling_factor,
    min_peak_height,
    list_intensity
  )
}

# Detect candidate peaks within an ROI, evaluate Gaussian peak shape and clean up overlapping peaks.
detect_evaluate_peaks_on_roi <- function(list_intensity_roi,
                                         rt_numbers_roi,
                                         min_peak_height,
                                         min_peak_ratio,
                                         min_fwhm,
                                         min_prominence_threshold,
                                         wlen,
                                         snr,
                                         peakshape,
                                         min_prominence_ratio,
                                         noise_level) {
  # min_peak_ratio and snr are preserved in the Python function signature, but are not used directly by this layer.
  invisible(min_peak_ratio)
  invisible(snr)
  list_peaks <- list()
  max_intensity <- max(list_intensity_roi)

  # The prominence threshold is simultaneously controlled by fixed threshold, ROI height, peak height threshold, and noise.
  roi_prominence <- max(
    min_prominence_threshold,
    max_intensity * min_prominence_ratio,
    0.33 * min_peak_height,
    2 * noise_level
  )

  # The first pass strictly uses wlen, and the parameter set is consistent with the Python scipy.find_peaks call.
  detected <- .find_peaks_scipy(
    list_intensity_roi,
    height = min_peak_height,
    distance = min_fwhm,
    prominence = roi_prominence,
    width = min_fwhm,
    wlen = wlen
  )

  # If there is no peak in the first round and the ROI is longer than wlen, cancel the prominence window and re-detect.
  if (length(detected$peaks) == 0L && length(list_intensity_roi) > wlen) {
    detected <- .find_peaks_scipy(
      list_intensity_roi,
      height = min_peak_height,
      distance = min_fwhm,
      prominence = roi_prominence,
      width = min_fwhm,
      wlen = NULL
    )
  }

  # The base span is checked peak-by-peak before calculating the Gaussian fit and the final JSON field.
  for (ii in seq_along(detected$peaks)) {
    if (detected$properties$right_bases[[ii]] -
        detected$properties$left_bases[[ii]] >= min_fwhm + 2L) {
      evaluated <- evaluate_roi_peak_json_(
        ii - 1L,
        list_intensity_roi,
        rt_numbers_roi,
        detected$peaks,
        detected$properties,
        peakshape,
        min_fwhm
      )
      if (!is.null(evaluated)) {
        list_peaks[[length(list_peaks) + 1L]] <- evaluated
      }
    }
  }

  check_overlap_peaks(list_peaks)
}

# Convert the ii-th ROI candidate peak into a named list with boundaries and peak shape score.
evaluate_roi_peak_json_ <- function(ii,
                                    list_intensity_roi,
                                    rt_numbers_roi,
                                    peaks,
                                    properties,
                                    peakshape,
                                    min_fwhm) {
  # ii comes from Python semantics, so accessing an R attribute vector adds 1 first.
  position <- as.integer(ii) + 1L
  left_index <- properties$left_bases[[position]]
  right_index <- properties$right_bases[[position]]

  # The apex and bounds of the Gaussian fit are relative to the current ROI and remain 0-based.
  gaussian <- evaluate_gaussian_peak_on_intensity_list(
    list_intensity_roi,
    properties$peak_heights[[position]],
    peaks[[position]],
    left_index,
    right_index
  )
  goodness_fitting <- gaussian[[1L]]
  sigma <- gaussian[[2L]]

  # This peak is retained only if the peak shape score strictly exceeds the threshold.
  if (goodness_fitting > peakshape) {
    halfwidth <- .python_int(max(abs(sigma), min_fwhm) * 6)
    proposed_left <- peaks[[position]] - halfwidth
    proposed_right <- peaks[[position]] + halfwidth
    left_index <- max(left_index, proposed_left)
    right_index <- min(right_index, proposed_right)

    # Maps the ROI internal position back to the scan number of the full mass track.
    left_base <- rt_numbers_roi[[left_index + 1L]]
    right_base <- rt_numbers_roi[[right_index + 1L]]
    apex <- rt_numbers_roi[[peaks[[position]] + 1L]]
    peak_area <- .python_int(sum(
      list_intensity_roi[seq.int(left_index + 1L, right_index + 1L)]
    ))

    return(list(
      apex = apex,
      peak_area = peak_area,
      height = .python_int(properties$peak_heights[[position]]),
      left_base = left_base,
      right_base = right_base,
      goodness_fitting = goodness_fitting
    ))
  }

  NULL
}

# Calculate cSelectivity for each peak on a mass track and scan points covered by all peaks.
`__peaks_cSelectivity_stats_` <- function(list_intensity, jpeaks) {
  # Summarize the closed interval between the left and right bases of each peak; unique corresponds to Python set.
  peak_datapoints <- unique(unlist(lapply(jpeaks, function(peak) {
    seq.int(peak$left_base, peak$right_base)
  }), use.names = FALSE))

  # Null peak lists must produce null intensities and null fractions without false subscripting.
  if (length(peak_datapoints) == 0L) {
    return(list(numeric(), integer()))
  }
  peak_intensities <- list_intensity[peak_datapoints + 1L]

  # For each peak, half the peak height is used as a threshold to compare the number of peak area points and the number of full track background points.
  selectivity <- vapply(jpeaks, function(peak) {
    threshold <- 0.5 * peak$height
    peak_level <- sum(peak_intensities > threshold)
    background_level <- sum(list_intensity > threshold)

    # Consistent with Python: scale is calculated only if 0 < peak_level <= background_level.
    if (background_level >= peak_level && peak_level > 0L) {
      peak_level / background_level
    } else {
      0
    }
  }, numeric(1))

  list(selectivity, sort(as.integer(peak_datapoints)))
}

# Standard Gaussian function, a is the height, mu is the center, and sigma is the standard deviation.
gaussian_function__ <- function(x, a, mu, sigma) {
  # The formula is exactly the same as Python a*np.exp(-(x-mu)**2/(2*sigma**2)).
  a * exp(-((x - mu)^2) / (2 * sigma^2))
}

# Use the coefficient of determination R^2 to evaluate the quality of the fit.
goodness_fitting__ <- function(y_orignal, y_fitted) {
  # Keep the original spelling of the Python parameter name y_oriignal to avoid interface drift.
  1 - sum((y_fitted - y_orignal)^2) /
    sum((y_orignal - mean(y_orignal))^2)
}

# Performs a Gaussian fit to the intensity points within the given bounds, returning R^2 and sigma.
evaluate_gaussian_peak_on_intensity_list <- function(intensity_list,
                                                     height,
                                                     apex,
                                                     left,
                                                     right) {
  # Inverse bounds produces empty array in Python, fit fails and initial sigma is NaN.
  if (right < left) return(list(0, NaN))

  # Python arange(left, right+1) includes the left and right boundaries, and the scan position remains 0-based.
  x <- seq.int(left, right)
  y <- intensity_list[x + 1L]
  initial_sigma <- sqrt(mean((x - mean(x))^2))
  goodness <- 0
  sigma <- initial_sigma

  # curve_fit will fail and return the fallback value of this function when the number of observation points is less than three parameters.
  if (length(x) < 3L) return(list(goodness, sigma))

  # Keep goodness=0 and initial sigma when fitting fails, consistent with the Python except branch.
  coefficients <- tryCatch(
    .fit_gaussian(x, y, c(height, apex, initial_sigma)),
    error = function(e) NULL
  )
  if (!is.null(coefficients)) {
    fitted <- gaussian_function__(
      x,
      coefficients[["a"]],
      coefficients[["mu"]],
      coefficients[["sigma"]]
    )
    goodness <- goodness_fitting__(y, fitted)
    sigma <- unname(coefficients[["sigma"]])
  }

  list(goodness, sigma)
}

# Use Gaussian fitting to estimate peak area within a given interval.
get_gaussian_peakarea_on_intensity_list <- function(intensity_list,
                                                    left,
                                                    right) {
  # Python will constrain right to length and use the half-open interval [left, right).
  right <- min(right, length(intensity_list))
  if (right <= left) return(0)
  x <- seq.int(left, right - 1L)
  y <- intensity_list[x + 1L]
  area <- 0

  # Use the highest point of the interval, the location mean and the overall standard deviation as the initial value of curve_fit.
  if (length(x) >= 3L) {
    initial_sigma <- sqrt(mean((x - mean(x))^2))
    coefficients <- tryCatch(
      .fit_gaussian(x, y, c(max(y), mean(x), initial_sigma)),
      error = function(e) NULL
    )

    # The total integral of the Gaussian is a * sqrt(2*pi*sigma^2).
    if (!is.null(coefficients)) {
      area <- coefficients[["a"]] * sqrt(
        2 * pi * coefficients[["sigma"]]^2
      )
    }
  }

  unname(area)
}

# Use LOWESS of chromatograms.R to smooth the entire mass track.
lowess_smooth_track <- function(list_intensity, number_of_scans) {
  fraction <- 0.05

  # Expand the local window when shorter than 200 scans, but keep the upper limit at 0.8.
  if (number_of_scans < 200) {
    fraction <- min(10 / number_of_scans, 0.8)
  }

  as.numeric(smooth_lowess(list_intensity, frac = fraction))
}

# Quickly detect if there happens to be a high and prominent elution peak.
quick_detect_unique_elution_peak <- function(
    intensity_track,
    min_peak_height = 100000,
    min_fwhm = 3,
    min_prominence_threshold_ratio = 0.2) {
  max_intensity <- max(intensity_track)
  prominence <- min_prominence_threshold_ratio * max_intensity
  unique_peak <- NULL

  # Python uses a strict greater than sign and does not enter detection when it is exactly equal to the lowest peak height.
  if (max_intensity > min_peak_height) {
    detected <- .find_peaks_scipy(
      intensity_track,
      height = min_peak_height,
      prominence = prominence,
      width = min_fwhm
    )

    # apex and height are returned only if one and only one peak is detected.
    if (length(detected$peaks) == 1L) {
      unique_peak <- list(
        apex = detected$peaks[[1L]],
        height = detected$properties$peak_heights[[1L]]
      )
    }
  }

  unique_peak
}

# Group peaks sorted by RT and clean up peaks that overlap with each other.
check_overlap_peaks <- function(list_peaks) {
  if (length(list_peaks) < 2L) return(list_peaks)

  # Continuous peaks are classified into the same cluster as long as they overlap with the last peak of the current cluster for more than 3 scans.
  clusters <- list()
  current_cluster <- list(list_peaks[[1L]])
  for (peak in list_peaks[-1L]) {
    if (`_check_overlap`(peak, tail(current_cluster, 1L)[[1L]])) {
      current_cluster[[length(current_cluster) + 1L]] <- peak
    } else {
      clusters[[length(clusters) + 1L]] <- current_cluster
      current_cluster <- list(peak)
    }
  }
  clusters[[length(clusters) + 1L]] <- current_cluster

  # Each cluster performs double peak boundaries splitting or multi-peak merging separately.
  unlist(lapply(clusters, cleanup_peak_cluster), recursive = FALSE)
}

# Determine whether the boundary overlap of two peaks strictly exceeds 3 scans.
`_check_overlap` <- function(peak1, peak2) {
  interval1 <- c(peak1$left_base, peak1$right_base)
  interval2 <- c(peak2$left_base, peak2$right_base)
  left_peak <- interval1
  right_peak <- interval2

  # Python distinguishes left peaks from right peaks by right_base size, rather than sorting by apex.
  if (interval1[[2L]] > interval2[[2L]]) {
    left_peak <- interval2
    right_peak <- interval1
  }
  overlap <- max(0, left_peak[[2L]] - right_peak[[1L]])
  overlap > 3
}

# Merge three or more overlapping peaks into one peak.
`_merge_peak_cluster` <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) return(cluster_peaks[[1L]])

  # Select the peak with the largest base span to inherit the other properties; when juxtaposed, the first one is retained.
  peak_sizes <- vapply(cluster_peaks, function(peak) {
    peak$right_base - peak$left_base
  }, numeric(1))
  largest <- cluster_peaks[[which.max(peak_sizes)]]

  # After merging, the boundaries are expanded to the leftmost and rightmost boundaries of the entire cluster.
  largest$left_base <- min(vapply(
    cluster_peaks, `[[`, numeric(1), "left_base"
  ))
  largest$right_base <- max(vapply(
    cluster_peaks, `[[`, numeric(1), "right_base"
  ))
  largest
}

# Correct the boundaries of an overlapping peak cluster.
cleanup_peak_cluster <- function(cluster_peaks) {
  if (length(cluster_peaks) == 1L) return(cluster_peaks)

  if (length(cluster_peaks) == 2L) {
    peak1 <- cluster_peaks[[1L]]
    peak2 <- cluster_peaks[[2L]]

    # When there are two peaks, sort the four bases and assign the first two and the last two to the two peaks respectively.
    bases <- sort(c(
      peak1$left_base, peak1$right_base,
      peak2$left_base, peak2$right_base
    ))
    peak1$left_base <- bases[[1L]]
    peak1$right_base <- bases[[2L]]
    peak2$left_base <- bases[[3L]]
    peak2$right_base <- bases[[4L]]

    # If the original apex is no longer inside the new boundary, use the boundary midpoint and truncate by int().
    if (!(peak1$left_base < peak1$apex && peak1$apex < peak1$right_base)) {
      peak1$apex <- .python_int(0.5 * (
        peak1$left_base + peak1$right_base
      ))
    }
    if (!(peak2$left_base < peak2$apex && peak2$apex < peak2$right_base)) {
      peak2$apex <- .python_int(0.5 * (
        peak2$left_base + peak2$right_base
      ))
    }

    return(list(peak1, peak2))
  }

  # Three or more overlapping peaks are merged into one peak according to Python rules.
  list(`_merge_peak_cluster`(cluster_peaks))
}

# Add up to three scan points to the left and right ends of the ROI that is too short.
extend_ROI <- function(ROI, number_of_scans) {
  # Only valid 0-based scan numbers greater than or equal to 0 are kept on the left.
  left <- c(ROI[[1L]] - 3L, ROI[[1L]] - 2L, ROI[[1L]] - 1L)
  left <- left[left >= 0L]

  # Only numbers strictly smaller than the total number of scans are kept on the right.
  last_scan <- tail(ROI, 1L)
  right <- c(last_scan + 1L, last_scan + 2L, last_scan + 3L)
  right <- right[right < number_of_scans]

  c(left, ROI, right)
}
