# Functions corresponding to asari/chromatograms.py.
#
# This file should contain mass track / EIC construction and RT alignment logic.
# The first implementation target is extract_mass_tracks().

# Read MS1 and MS2 scans from mzML files.
# MS1 data are intensity filtered, m/z binning, and merged to generate mass tracks; MS2 spectra are saved separately.
extract_mass_tracks <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  rt_times <- numeric()
  ms2_spectra <- list()
  ii <- 0L

  ms_expt <- mzR::openMSfile(infile)
  on.exit(mzR::close(ms_expt), add = TRUE)
  scan_header <- mzR::header(ms_expt)

  # Calling peaks() on mzR scan by scan will repeatedly locate and decode the same mzML, which will cause huge overhead on real files.
  # Read all scans at once and then access them in the original order. The returned m/z and intensity matrices are the same as reading them sequentially.
  scan_peaks <- mzR::peaks(ms_expt)
  # First press scan to save the data points, and finally flatten them all at once; avoid expanding a large list every time a point is added.
  scan_points <- vector("list", nrow(scan_header))

  timestamp <- tryCatch(
    as.numeric(as.POSIXct(
      mzR::runInfo(ms_expt)$startTimeStamp,
      format = "%Y-%m-%dT%H:%M:%SZ",
      tz = "UTC"
    )),
    error = function(e) NULL
  )

  for (spec in seq_len(nrow(scan_header))) {
    if (scan_header$msLevel[[spec]] == 1L) {
      rt_times <- c(rt_times, scan_header$retentionTime[[spec]])
      spectrum <- scan_peaks[[spec]]
      intensities <- trunc(spectrum[, "intensity"])
      good_positions <- intensities > min_intensity

      # Construct all triples of the current scan at once, maintaining the original m/z order and 0-based scan number.
      positions <- which(good_positions)
      if (length(positions) > 0L) {
        scan_points[[spec]] <- Map(
          function(mz, intensity) c(mz, ii, intensity),
          spectrum[positions, "mz"],
          intensities[positions]
        )
      }
      ii <- ii + 1L
    } else if (scan_header$msLevel[[spec]] == 2L) {
      precursor_mz <- scan_header$precursorMZ[[spec]]
      if (is.na(precursor_mz)) next

      spectrum <- scan_peaks[[spec]]
      keep <- spectrum[, "intensity"] >= min_intensity &
        spectrum[, "mz"] < precursor_mz - 1
      peaks <- lapply(which(keep), function(jj) {
        c(spectrum[jj, "mz"], spectrum[jj, "intensity"])
      })

      if (length(peaks) > 0L) {
        ms2_spectra[[length(ms2_spectra) + 1L]] <- list(
          precursor_mz = precursor_mz,
          rtime = scan_header$retentionTime[[spec]],
          peaks = peaks
        )
      }
    }
  }

  # Flatten and group by thousandth m/z at one time to avoid duplicating the list when repeatedly appending to each bin.
  alldata <- unlist(scan_points, recursive = FALSE, use.names = FALSE)
  mz_keys <- vapply(
    alldata,
    function(point) as.integer(point[[1L]] * 1000),
    integer(1)
  )
  mz_tree <- split(
    alldata,
    factor(mz_keys, levels = unique(mz_keys)),
    drop = TRUE
  )

  rm(alldata, scan_points)
  rt_numbers <- seq_len(length(rt_times)) - 1L
  rt_length <- length(rt_numbers)

  good_bins <- get_thousandth_bins(
    mz_tree, mz_tolerance_ppm, min_timepoints, min_peak_height
  )
  tracks <- list()
  for (bin in good_bins) {
    tracks <- c(tracks, bin_to_mass_tracks(bin, rt_length, mz_tolerance_ppm))
  }

  tracks <- tracks[order(vapply(tracks, function(x) x[[1L]], numeric(1)))]
  merged <- list()
  to_remove <- integer()
  tracks_to_merge <- check_close_mzs(
    vapply(tracks, function(x) x[[1L]], numeric(1)),
    mz_tolerance_ppm
  )
  for (pair in tracks_to_merge) {
    a <- pair[[1L]]
    b <- pair[[2L]]
    merged[[length(merged) + 1L]] <- merge_two_mass_tracks(
      tracks[[a]], tracks[[b]]
    )
    to_remove <- c(to_remove, a, b)
  }
  updated_tracks <- c(
    tracks[!seq_along(tracks) %in% to_remove],
    merged
  )

  list(
    rt_numbers = rt_numbers,
    rt_times = rt_times,
    tracks = updated_tracks,
    ms2_spectra = ms2_spectra,
    acquisition_time = timestamp
  )
}

# Corresponds to the original name of Python extract_massTracks_; retains its uppercase and lowercase and trailing underscore interfaces.
extract_massTracks_ <- function(infile,
                                mz_tolerance_ppm = 5,
                                min_intensity = 100,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  extract_mass_tracks(
    infile, mz_tolerance_ppm, min_intensity, min_timepoints, min_peak_height
  )
}

# Arrange data points within a limited m/z range into a mass track covering the full RT range.
# Represents m/z calculated combining the median and highest intensity points; when multiple points are present on the same scan, the highest intensity is retained.
extract_single_track_full_rt_length <- function(bin,
                                                rt_length,
                                                intensity_data_type = numeric) {
  mzs <- vapply(bin, function(x) x[[1L]], numeric(1))
  ints <- vapply(bin, function(x) x[[3L]], numeric(1))
  mz <- 0.5 * (stats::median(mzs) + mzs[[which.max(ints)]])
  intensity_track <- intensity_data_type(rt_length)

  for (r in bin) {
    scan_index <- r[[2L]] + 1L
    intensity_track[[scan_index]] <- max(
      r[[3L]], intensity_track[[scan_index]]
    )
  }

  list(mz, intensity_track)
}

# Corresponds to the original name of Python extract_single_track_fullrt_length.
extract_single_track_fullrt_length <- function(
    bin, rt_length, INTENSITY_DATA_TYPE = numeric) {
  extract_single_track_full_rt_length(bin, rt_length, INTENSITY_DATA_TYPE)
}

# Determine whether an m/z bin should constitute one or multiple mass tracks.
# When the m/z range is small, one line is constructed directly; when the range is large, m/z is first clustered and then constructed group by group.
bin_to_mass_tracks <- function(bin_data_tuples, rt_length, mz_tolerance_ppm = 5) {
  bin_data_tuples <- bin_data_tuples[
    order(vapply(bin_data_tuples, function(x) x[[1L]], numeric(1)))
  ]
  mz_range <- bin_data_tuples[[length(bin_data_tuples)]][[1L]] -
    bin_data_tuples[[1L]][[1L]]
  mz_tolerance <- bin_data_tuples[[1L]][[1L]] * 0.000001 * mz_tolerance_ppm

  if (mz_range < mz_tolerance * 2) {
    return(list(
      extract_single_track_full_rt_length(bin_data_tuples, rt_length)
    ))
  }

  rois <- build_chromatogram_by_mz_clustering(
    bin_data_tuples, rt_length, mz_tolerance
  )
  lapply(rois, function(roi) {
    extract_single_track_full_rt_length(roi, rt_length)
  })
}

# Corresponds to build_chromatogram_intensity_aware: starting from the highest intensity seed, absorbing points within the m/z tolerance group by group.
build_chromatogram_intensity_aware <- function(
    bin_data_tuples, rt_length, mz_tolerance) {
  invisible(rt_length)
  if (length(bin_data_tuples) == 0L) return(list())
  intensity <- vapply(bin_data_tuples, function(point) point[[3L]], numeric(1))
  remaining <- bin_data_tuples[order(intensity, decreasing = TRUE)]
  assigned <- list()
  while (length(remaining)) {
    seed <- list(remaining[[1L]])
    next_remaining <- list()
    if (length(remaining) > 1L) {
      for (point in remaining[-1L]) {
        if (abs(point[[1L]] - seed[[1L]][[1L]]) < mz_tolerance) {
          seed[[length(seed) + 1L]] <- point
        } else {
          next_remaining[[length(next_remaining) + 1L]] <- point
        }
      }
    }
    assigned[[length(assigned) + 1L]] <- seed
    remaining <- next_remaining
  }
  assigned
}

# Give the m/z data that needs to be split to the nearest neighbor clustering function based on m/z seeds.
# rt_length is currently not involved in calculations and is only reserved to correspond to Python function signatures.
build_chromatogram_by_mz_clustering <- function(bin_data_tuples,
                                                rt_length,
                                                mz_tolerance) {
  nn_cluster_by_mz_seeds(
    bin_data_tuples,
    mz_tolerance,
    presorted = FALSE
  )
}

# Merge two mass tracks belonging to the same ion signal.
# Representative m/z of the two tracks were averaged, and the intensities over the full RT range were summed position by position.
merge_two_mass_tracks <- function(t1, t2) {
  list(
    0.5 * (t1[[1L]] + t2[[1L]]),
    t1[[2L]] + t2[[2L]]
  )
}

# Corresponds to get_thousandth_bins internal __rough_check_consecutive_scans__.
`__rough_check_consecutive_scans__` <- function(
    datatuples, gap_allowed = 2, min_timepoints = 5) {
  checked <- TRUE
  check_max_len <- 4 * min_timepoints
  if (length(datatuples) < check_max_len) {
    min_check_val <- gap_allowed + min_timepoints - 1
    rts <- sort(vapply(datatuples, function(point) point[[2L]], numeric(1)))
    stops <- seq.int(min_timepoints, length(rts))
    steps <- vapply(
      stops,
      function(index) rts[[index]] - rts[[index - min_timepoints + 1L]],
      numeric(1)
    )
    if (min(steps) > min_check_val) checked <- FALSE
  }
  checked
}

# Corresponds to get_thousandth_bins internal __check_min_peak_height__.
`__check_min_peak_height__` <- function(datatuples, min_peak_height) {
  max(vapply(datatuples, function(point) point[[3L]], numeric(1))) >= min_peak_height
}

# Merge adjacent or ppm-range thousandth m/z bins and filter weak and scattered signals.
# The returned good bins will continue to be handed over to bin_to_mass_tracks() to build mass tracks.
get_thousandth_bins <- function(mz_tree,
                                mz_tolerance_ppm = 5,
                                min_timepoints = 5,
                                min_peak_height = 1000) {
  tolerance <- 0.000001 * mz_tolerance_ppm
  keys <- sort(as.integer(
    names(mz_tree)[lengths(mz_tree) >= min_timepoints]
  ))
  bins_of_bins <- list()
  current_bin <- keys[[1L]]

  if (length(keys) > 1L) {
    for (ii in seq.int(2L, length(keys))) {
      delta <- keys[[ii]] - keys[[ii - 1L]]
      if (delta == 1L || delta < tolerance * keys[[ii]]) {
        current_bin <- c(current_bin, keys[[ii]])
      } else {
        bins_of_bins[[length(bins_of_bins) + 1L]] <- current_bin
        current_bin <- keys[[ii]]
      }
    }
  }

  bins_of_bins[[length(bins_of_bins) + 1L]] <- current_bin
  good_bins <- list()

  for (bin in bins_of_bins) {
    datatuples <- list()
    for (key in bin) {
      datatuples <- c(datatuples, mz_tree[[as.character(key)]])
    }

    if (`__check_min_peak_height__`(datatuples, min_peak_height) &&
        `__rough_check_consecutive_scans__`(
          datatuples, min_timepoints = min_timepoints
        )) {
      good_bins[[length(good_bins) + 1L]] <- datatuples
    }
  }

  good_bins
}

# Corresponds to linear interpolation and extrapolation of SciPy interp1d(..., fill_value = "extrapolate").
linear_interpolate_with_extrapolation <- function(x, y, xout) {
  interpolation_order <- order(x)
  x <- x[interpolation_order]
  y <- y[interpolation_order]
  result <- stats::approx(x, y, xout = xout)$y

  left <- xout < x[[1L]]
  result[left] <- y[[1L]] + (xout[left] - x[[1L]]) *
    (y[[2L]] - y[[1L]]) / (x[[2L]] - x[[1L]])

  last <- length(x)
  right <- xout > x[[last]]
  result[right] <- y[[last]] + (xout[right] - x[[last]]) *
    (y[[last]] - y[[last - 1L]]) / (x[[last]] - x[[last - 1L]])

  result
}

# Use the matching landmark peaks to fit LOWESS and map the sample scan to the reference sample scan.
# Returns both forward and reverse mappings, saving only scans that have changed and have not exceeded the true range.
rt_lowess_calibration <- function(good_landmark_peaks,
                                  selected_reference_landmark_peaks,
                                  sample_rt_numbers,
                                  reference_rt_numbers,
                                  num_iterations,
                                  sample_name,
                                  outdir) {
  reference_rt_bound <- max(reference_rt_numbers)
  sample_rt_bound <- max(sample_rt_numbers)
  rt_rightend <- 1.1 * sample_rt_bound
  xx <- rep(-0.1 * sample_rt_bound, 3L)
  yy <- rep(-0.1 * sample_rt_bound, 3L)

  rt_cal <- clean_rt_calibration_points(Map(
    function(sample_peak, reference_peak) {
      c(sample_peak$apex, reference_peak$apex)
    },
    good_landmark_peaks,
    selected_reference_landmark_peaks
  ))

  xx <- c(
    xx,
    vapply(rt_cal, function(pair) pair[[1L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  yy <- c(
    yy,
    vapply(rt_cal, function(pair) pair[[2L]], numeric(1)),
    rep(rt_rightend, 3L)
  )

  frac <- 0.6 - 0.004 * (length(rt_cal) - 50)
  frac <- max(0.2, min(frac, 0.6))

  lowess_predicted <- hacked_lowess(
    yy,
    xx,
    frac = frac,
    it = num_iterations,
    xvals = sample_rt_numbers
  )

  ref_interpolated <- linear_interpolate_with_extrapolation(
    lowess_predicted,
    sample_rt_numbers,
    reference_rt_numbers
  )

  lowess_predicted <- as.integer(round(lowess_predicted))
  keep_forward <- sample_rt_numbers != lowess_predicted &
    lowess_predicted >= 0L &
    lowess_predicted <= reference_rt_bound
  rt_cal_dict <- stats::setNames(
    lowess_predicted[keep_forward],
    sample_rt_numbers[keep_forward]
  )

  ref_interpolated <- as.integer(round(ref_interpolated))
  keep_reverse <- reference_rt_numbers != ref_interpolated &
    ref_interpolated >= 0L &
    ref_interpolated <= sample_rt_bound
  reverse_rt_cal_dict <- stats::setNames(
    ref_interpolated[keep_reverse],
    reference_rt_numbers[keep_reverse]
  )

  list(rt_cal_dict, reverse_rt_cal_dict)
}

# Duplicate RT landmark pairs were removed and outliers with deviations exceeding three times the population standard deviation were filtered.
clean_rt_calibration_points <- function(rt_cal_pairs) {
  pair_keys <- vapply(
    rt_cal_pairs,
    function(pair) paste(pair, collapse = "\r"),
    character(1)
  )
  unique_pairs <- rt_cal_pairs[!duplicated(pair_keys)]
  deltas <- vapply(
    unique_pairs,
    function(pair) pair[[1L]] - pair[[2L]],
    numeric(1)
  )
  mean_delta <- mean(deltas)
  std3x <- sqrt(mean((deltas - mean_delta)^2)) * 3
  low <- mean_delta - std3x
  high <- mean_delta + std3x

  kept <- unique_pairs[deltas > low & deltas < high]
  kept[order(
    vapply(kept, function(pair) pair[[1L]], numeric(1)),
    vapply(kept, function(pair) pair[[2L]], numeric(1))
  )]
}

# Fit using R's LOWESS and linearly interpolate the fitted values to the specified scan coordinates.
hacked_lowess <- function(yy, xx, frac, it, xvals) {
  fitted <- stats::lowess(
    xx,
    yy,
    f = frac,
    iter = it,
    delta = 0
  )
  stats::approx(
    fitted$x,
    fitted$y,
    xout = xvals,
    ties = "ordered"
  )$y
}

# Corresponds to Python’s original name __hacked_lowess__.
`__hacked_lowess__` <- function(yy, xx, frac, it, xvals) {
  hacked_lowess(yy, xx, frac, it, xvals)
}

# Debug version of RT calibration: returns the same mapping and outputs landmark and calibration curve plots.
rt_lowess_calibration_debug <- function(good_landmark_peaks,
                                         selected_reference_landmark_peaks,
                                         sample_rt_numbers,
                                         reference_rt_numbers,
                                         num_iterations,
                                         sample_name,
                                         outdir) {
  calibration <- rt_lowess_calibration(
    good_landmark_peaks,
    selected_reference_landmark_peaks,
    sample_rt_numbers,
    reference_rt_numbers,
    num_iterations,
    sample_name,
    outdir
  )

  sample_rt_bound <- max(sample_rt_numbers)
  rt_rightend <- 1.1 * sample_rt_bound
  rt_cal <- clean_rt_calibration_points(Map(
    function(sample_peak, reference_peak) {
      c(sample_peak$apex, reference_peak$apex)
    },
    good_landmark_peaks,
    selected_reference_landmark_peaks
  ))
  xx <- c(
    rep(-0.1 * sample_rt_bound, 3L),
    vapply(rt_cal, function(pair) pair[[1L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  yy <- c(
    rep(-0.1 * sample_rt_bound, 3L),
    vapply(rt_cal, function(pair) pair[[2L]], numeric(1)),
    rep(rt_rightend, 3L)
  )
  frac <- 0.6 - 0.004 * (length(rt_cal) - 50)
  frac <- max(0.2, min(frac, 0.6))
  lowess_predicted <- as.integer(round(hacked_lowess(
    yy,
    xx,
    frac = frac,
    it = num_iterations,
    xvals = sample_rt_numbers
  )))

  plot_file <- file.path(
    outdir,
    "export",
    paste0(sample_name, "_rtime_alignment_result.png")
  )
  grDevices::png(plot_file)
  graphics::plot(
    xx[4L:(length(xx) - 3L)],
    yy[4L:(length(yy) - 3L)],
    xlab = "Sample Retention Time (sec)",
    ylab = "Reference Retention Time (sec)",
    main = "Retention Time Alignment",
    col = "black",
    pch = 1
  )
  graphics::lines(sample_rt_numbers, lowess_predicted, col = "red")
  graphics::legend(
    "topleft",
    legend = c("Landmark peaks", "Alignment function"),
    col = c("black", "red"),
    pch = c(1, NA),
    lty = c(NA, 1)
  )
  grDevices::dev.off()

  calibration
}

# Savitzky-Golay spline placeholder function not yet implemented in Python.
savitzky_golay_spline <- function(good_landmark_peaks,
                                   selected_reference_landmark_peaks,
                                   sample_rt_numbers,
                                   reference_rt_numbers) {
  invisible(NULL)
}

# DWT RT calibration placeholder function in Python explicitly marked as not yet implemented.
dwt_rt_calibrate <- function(good_landmark_peaks,
                             selected_reference_landmark_peaks,
                             sample_rt_numbers,
                             reference_rt_numbers) {
  stop("Not implemented", call. = FALSE)
}

# Copy the intensity value to the calibrated scan position according to the RT calibration map.
remap_intensity_track <- function(intensity_track, new, rt_cal_dict) {
  new[seq_along(intensity_track)] <- intensity_track

  for (ii in seq_along(rt_cal_dict)) {
    source_index <- as.integer(names(rt_cal_dict)[[ii]]) + 1L
    target_index <- as.integer(round(rt_cal_dict[[ii]])) + 1L
    if (source_index >= 1L && source_index <= length(intensity_track) &&
        target_index >= 1L && target_index <= length(new)) {
      new[[target_index]] <- intensity_track[[source_index]]
    }
  }

  new
}

# Fill the window with the nearest boundary value and do a simple moving average of the intensity track.
smooth_moving_average <- function(list_intensity, size = 9) {
  left_width <- floor(size / 2)
  right_width <- size - left_width - 1L

  vapply(seq_along(list_intensity), function(ii) {
    positions <- seq.int(ii - left_width, ii + right_width)
    positions <- pmin(length(list_intensity), pmax(1L, positions))
    mean(list_intensity[positions])
  }, numeric(1))
}

# Use LOWESS to smooth noisy intensity tracks.
smooth_lowess <- function(list_intensity, frac = 0.02) {
  scan_numbers <- seq_along(list_intensity) - 1L
  hacked_lowess(
    list_intensity,
    scan_numbers,
    frac = frac,
    it = 1,
    xvals = scan_numbers
  )
}
