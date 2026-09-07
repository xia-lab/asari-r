# Corresponds to Python asari/constructors.py: Complete R implementation of MassGrid and CompositeMap.
#
# Data structure convention:
# 1. Python mutable objects use environment in R to ensure that the internal modifications of the method can be retained.
# 2. MassGrid uses data.frame; the first column is mz, and the remaining columns save 0-based track IDs according to sample names.
# 3. Python’s external scan, feature and track IDs continue to remain 0-based; 1 is added when accessing the R vector.
# 4. Python None / pandas NaN is represented by NA in MassGrid.

# Reads the fields of an environment or list, allowing default values to be provided for missing fields.
.constructors_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && !is.null(object[[name]])) return(object[[name]])
  default
}

# Modify mutable object fields; Python objects that need to be modified in constructors must use environment.
.constructors_set <- function(object, name, value) {
  if (!is.environment(object)) {
    stop("Mutable constructor objects must be environments.", call. = FALSE)
  }
  assign(name, value, envir = object)
  invisible(value)
}

# Shallow copy a common named list into an environment to facilitate testing and subsequent access to samples.R.
.constructors_as_environment <- function(object, class_name = NULL) {
  if (is.environment(object)) {
    if (!is.null(class_name) && !inherits(object, class_name)) {
      class(object) <- unique(c(class_name, class(object)))
    }
    return(object)
  }
  if (!is.list(object)) stop("object must be a list or environment.", call. = FALSE)
  result <- list2env(object, parent = emptyenv())
  if (!is.null(class_name)) class(result) <- c(class_name, "environment")
  result
}

# Read sample registry entries according to Python registry[id] semantics, compatible with named IDs and 0-based numeric IDs.
.constructors_registry_get <- function(registry, sample_id) {
  if (is.environment(registry)) {
    key <- as.character(sample_id)
    if (exists(key, envir = registry, inherits = FALSE)) {
      return(get(key, envir = registry, inherits = FALSE))
    }
  }
  if (!is.list(registry)) stop("sample_registry must be a list or environment.", call. = FALSE)
  key <- as.character(sample_id)
  if (!is.null(names(registry)) && key %in% names(registry)) return(registry[[key]])
  if (is.numeric(sample_id)) return(registry[[as.integer(sample_id) + 1L]])
  stop("Unknown sample id: ", key, call. = FALSE)
}

# When samples.py has not yet been migrated, create a SimpleSample from the sample_factory or registry provided by the experiment.
.constructors_make_simple_sample <- function(registry,
                                             experiment,
                                             database_mode,
                                             mode,
                                             is_reference = FALSE) {
  factory <- .constructors_get(experiment, "sample_factory")
  if (is.function(factory)) {
    return(factory(
      registry,
      experiment = experiment,
      database_mode = database_mode,
      mode = mode,
      is_reference = is_reference
    ))
  }

  # If the registration item itself is already a SimpleSample environment, directly add the reference and status.
  sample <- .constructors_as_environment(registry, "SimpleSample")
  sample$experiment <- experiment
  sample$database_mode <- database_mode
  sample$mode <- mode
  sample$is_reference <- is_reference
  if (is.null(sample$is_rt_aligned)) sample$is_rt_aligned <- is_reference
  if (is.null(sample$rt_cal_dict)) sample$rt_cal_dict <- numeric()
  if (is.null(sample$reverse_rt_cal_dict)) sample$reverse_rt_cal_dict <- numeric()
  if (is.null(sample$list_scan_numbers) && !is.null(sample$rt_numbers)) {
    sample$list_scan_numbers <- sample$rt_numbers
  }
  sample
}

# Corresponds to SimpleSample.get_masstracks_and_anchors / get_mass_tracks_for_sample.
.constructors_get_mass_tracks <- function(sample) {
  method <- .constructors_get(sample, "get_masstracks_and_anchors")
  if (is.function(method)) return(method())
  method <- .constructors_get(sample, "get_mass_tracks_for_sample")
  if (is.function(method)) return(method(sample))
  tracks <- .constructors_get(sample, "list_mass_tracks")
  if (!is.null(tracks)) return(tracks)
  tracks <- .constructors_get(sample, "mass_tracks")
  if (!is.null(tracks)) return(tracks)
  stop("Sample does not provide mass tracks.", call. = FALSE)
}

# Copy linear interpolation from scipy.interpolate.interp1d(..., fill_value="extrapolate").
.constructors_interp_extrapolate <- function(x, y, xout) {
  if (length(x) < 2L || length(y) != length(x)) {
    stop("Interpolation requires at least two paired points.", call. = FALSE)
  }
  if (exists("linear_interpolate_with_extrapolation", mode = "function")) {
    return(linear_interpolate_with_extrapolation(x, y, xout))
  }

  order_x <- order(x)
  x <- as.numeric(x[order_x])
  y <- as.numeric(y[order_x])
  result <- stats::approx(x, y, xout = xout, ties = "ordered")$y
  left <- xout < x[[1L]]
  right <- xout > x[[length(x)]]
  result[left] <- y[[1L]] + (xout[left] - x[[1L]]) *
    (y[[2L]] - y[[1L]]) / (x[[2L]] - x[[1L]])
  last <- length(x)
  result[right] <- y[[last]] + (xout[right] - x[[last]]) *
    (y[[last]] - y[[last - 1L]]) / (x[[last]] - x[[last - 1L]])
  result
}

# Read an R vector of scan-numbered names, simulating Python dict.get(key, default).
.constructors_dict_get <- function(dictionary, key, default = NULL) {
  if (is.null(dictionary) || length(dictionary) == 0L) return(default)
  key <- as.character(key)
  if (!is.null(names(dictionary)) && key %in% names(dictionary)) {
    return(unname(dictionary[[key]]))
  }
  default
}

# Write the values in the named map to the vector names=key, retaining the last assignment.
.constructors_dict_set <- function(dictionary, key, value) {
  if (is.null(dictionary)) dictionary <- numeric()
  dictionary[[as.character(key)]] <- value
  dictionary
}

# Reproduce scipy.ndimage.maximum_filter1d(size=2, mode="constant").
.constructors_maximum_filter_size2 <- function(values) {
  values <- as.numeric(values)
  if (length(values) == 0L) return(values)
  pmax(c(0, values[-length(values)]), values)
}

# Convert Python [start:stop] half-open slices to R subscripts, and handle negative bounds and out-of-bounds truncation.
.constructors_python_slice <- function(values, start, stop) {
  n <- length(values)
  start <- as.integer(start)
  stop <- as.integer(stop)
  if (start < 0L) start <- max(n + start, 0L)
  if (stop < 0L) stop <- max(n + stop, 0L)
  start <- min(max(start, 0L), n)
  stop <- min(max(stop, 0L), n)
  if (stop <= start) return(values[0])
  values[seq.int(start + 1L, stop)]
}

# Corresponds to mass2chem.find_mzdiff_pairs_from_masstracks, finding 13C/12C and Na/H anchors.
.constructors_find_mzdiff_pairs <- function(list_mass_tracks,
                                            list_mz_diff = c(1.003355, 21.9820),
                                            mz_tolerance_ppm = 5) {
  if (length(list_mass_tracks) == 0L) return(list())
  mzs <- vapply(list_mass_tracks, `[[`, numeric(1), "mz")
  pairs <- list()

  # Find the mass track with the smallest absolute error within the tolerance for each target difference.
  for (difference in list_mz_diff) {
    for (track in list_mass_tracks) {
      target <- track$mz + difference
      errors <- abs(mzs - target)
      best <- which.min(errors)
      if (errors[[best]] < target * mz_tolerance_ppm * 1e-6) {
        pairs[[length(pairs) + 1L]] <- c(
          track$id_number,
          list_mass_tracks[[best]]$id_number
        )
      }
    }
  }
  pairs
}

# Get Python MassGrid[row, sample]: both row_id and track ID remain 0-based.
.constructors_grid_get <- function(grid, sample_name, row_id) {
  if (!(sample_name %in% names(grid))) return(NA_real_)
  row <- as.integer(row_id) + 1L
  if (row < 1L || row > nrow(grid)) return(NA_real_)
  grid[[sample_name]][[row]]
}

# -----------------------------------------------------------------------------
# Python class MassGrid: 7 methods and 1 nested helper function
# -----------------------------------------------------------------------------

# Corresponds to MassGrid.__init__: Creates a variable MassGrid environment.
MassGrid__init__ <- function(cmap, experiment) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("MassGrid", "environment")
  self$experiment <- experiment
  self$CMAP <- cmap
  self$reference_sample_instance <- cmap$reference_sample_instance
  self$max_ref_rtime <- cmap$max_ref_rtime
  self$list_sample_names <- cmap$list_sample_names
  self$`_number_of_samples_` <- cmap$`_number_of_samples_`
  self$MassGrid <- NULL
  self$`_mz_landmarks_` <- integer()
  self$anchor_mz_pairs <- list()
  self
}

# Corresponds to MassGrid.build_grid_sample_wise: adds dynamic reference m/z grid sample by sample.
MassGrid_build_grid_sample_wise <- function(self) {
  MassGrid__initiate_mass_grid(self)
  sample_ids <- self$experiment$valid_sample_ids
  sample_ids <- sample_ids[sample_ids != self$experiment$reference_sample_id]

  # Python creates a SimpleSample for each non-reference entry and calls add_sample.
  for (sample_id in sample_ids) {
    sample <- .constructors_make_simple_sample(
      .constructors_registry_get(self$experiment$sample_registry, sample_id),
      self$experiment,
      self$experiment$database_mode,
      self$experiment$mode,
      is_reference = FALSE
    )
    MassGrid_add_sample(self, sample)
  }
  invisible(self)
}

# Corresponds to MassGrid.build_grid_by_centroiding: cluster the track m/z of all samples at once.
MassGrid_build_grid_by_centroiding <- function(self) {
  all_tracks <- list()

  # The tuple structure holds (mz, 0-based track_id, 0-based sample_index).
  for (sample_index in seq_len(self$`_number_of_samples_`) - 1L) {
    sample_id <- self$experiment$valid_sample_ids[[sample_index + 1L]]
    registry <- .constructors_registry_get(
      self$experiment$sample_registry, sample_id
    )
    track_mzs <- .constructors_get(registry, "track_mzs", list())
    for (track_info in track_mzs) {
      all_tracks[[length(all_tracks) + 1L]] <- c(
        track_info[[1L]], track_info[[2L]], sample_index
      )
    }
  }

  # Python tuple sort sorts by mz, track_id, sample_index.
  ordering <- order(
    vapply(all_tracks, `[[`, numeric(1), 1L),
    vapply(all_tracks, `[[`, numeric(1), 2L),
    vapply(all_tracks, `[[`, numeric(1), 3L)
  )
  all_tracks <- all_tracks[ordering]
  all_bins <- MassGrid_bin_track_mzs(
    self, all_tracks, self$experiment$reference_sample_id
  )

  # Initialize the sample column and fill in the track ID in each bin into the corresponding sample column.
  grid <- as.data.frame(
    matrix(NA_real_, nrow = length(all_bins), ncol = self$`_number_of_samples_`),
    stringsAsFactors = FALSE
  )
  names(grid) <- self$CMAP$list_sample_names
  for (bin_index in seq_along(all_bins)) {
    for (track_tuple in all_bins[[bin_index]][[2L]]) {
      sample_column <- as.integer(track_tuple[[3L]]) + 1L
      grid[[sample_column]][[bin_index]] <- track_tuple[[2L]]
    }
  }
  mz_list <- vapply(all_bins, `[[`, numeric(1), 1L)
  self$MassGrid <- cbind(mz = mz_list, grid)
  names(self$MassGrid) <- c("mz", self$CMAP$list_sample_names)

  # Calculate m/z difference anchor points and convert 0-based pair into R 1-based landmark row positions.
  list_mass_tracks <- lapply(seq_along(mz_list), function(ii) {
    list(id_number = ii - 1L, mz = mz_list[[ii]])
  })
  self$anchor_mz_pairs <- .constructors_find_mzdiff_pairs(
    list_mass_tracks,
    mz_tolerance_ppm = self$experiment$parameters$mz_tolerance_ppm
  )
  flattened <- flatten_tuplelist(self$anchor_mz_pairs)
  self$`_mz_landmarks_` <- if (is.null(flattened)) integer() else {
    sort(as.integer(flattened) + 1L)
  }

  # Create reference and remaining SimpleSample instances, keeping experiment.all_samples order.
  self$reference_sample_instance$rt_cal_dict <- numeric()
  self$reference_sample_instance$reverse_rt_cal_dict <- numeric()
  self$experiment$all_samples[[length(self$experiment$all_samples) + 1L]] <-
    self$reference_sample_instance
  for (sample_id in self$experiment$valid_sample_ids) {
    if (sample_id == self$experiment$reference_sample_id) next
    sample <- .constructors_make_simple_sample(
      .constructors_registry_get(self$experiment$sample_registry, sample_id),
      self$experiment,
      self$experiment$database_mode,
      self$experiment$mode,
      is_reference = FALSE
    )
    self$experiment$all_samples[[length(self$experiment$all_samples) + 1L]] <- sample
  }
  invisible(self)
}

# Corresponds to MassGrid._initiate_mass_grid: Use reference sample to create the first version of MassGrid.
MassGrid__initiate_mass_grid <- function(self) {
  reference_sample <- self$reference_sample_instance
  reference_sample$rt_cal_dict <- numeric()
  reference_sample$reverse_rt_cal_dict <- numeric()
  reference_tracks <- reference_sample$list_mass_tracks
  self$`_mz_landmarks_` <- sort(as.integer(reference_sample$`_mz_landmarks_`))
  reference_mzs <- vapply(reference_tracks, `[[`, numeric(1), "mz")

  # All sample columns are filled with NA first, and the reference sample column is written with Python 0-based id_number.
  grid <- data.frame(mz = reference_mzs, check.names = FALSE)
  for (sample_name in self$list_sample_names) grid[[sample_name]] <- NA_real_
  grid[[reference_sample$name]] <- vapply(
    reference_tracks, `[[`, numeric(1), "id_number"
  )
  self$MassGrid <- grid
  self$experiment$all_samples[[length(self$experiment$all_samples) + 1L]] <-
    reference_sample
  invisible(self)
}

# Corresponds to MassGrid.add_sample: add a new sample using landmark guidance.
MassGrid_add_sample <- function(self, sample) {
  message("Adding sample to MassGrid, ", sample$name)
  mz_list <- vapply(sample$track_mzs, `[[`, numeric(1), 1L)

  # The landmark interface and mapping position of mass_functions.R use R 1-based.
  mapping <- landmark_guided_mapping(
    self$MassGrid$mz,
    as.integer(self$`_mz_landmarks_`),
    mz_list,
    as.integer(sample$`_mz_landmarks_`),
    std_ppm = self$experiment$parameters$mz_tolerance_ppm,
    correction_tolerance_ppm =
      self$experiment$parameters$correction_tolerance_ppm
  )

  # The new grid copies the old row; the mapping sample position is decremented by 1 and saved as the Python track ID.
  new_grid <- data.frame(mz = mapping$new_reference_mzlist, check.names = FALSE)
  for (sample_name in self$list_sample_names) new_grid[[sample_name]] <- NA_real_
  old_rows <- seq_len(nrow(self$MassGrid))
  for (column in names(self$MassGrid)[-1L]) {
    new_grid[[column]][old_rows] <- self$MassGrid[[column]]
  }
  new_map <- mapping$new_reference_map2
  new_grid[[sample$name]] <- ifelse(is.na(new_map), NA_real_, new_map - 1L)
  self$MassGrid <- new_grid
  self$`_mz_landmarks_` <- as.integer(mapping$REF_landmarks)
  sample$mz_calibration_ratio <- mapping$correction_ratio
  self$experiment$all_samples[[length(self$experiment$all_samples) + 1L]] <- sample
  invisible(self)
}

# Corresponds to MassGrid.bin_track_mzs inline def __get_bin__.
MassGrid_bin_track_mzs__get_bin <- function(bin_data_tuples) {
  list(
    stats::median(vapply(bin_data_tuples, `[[`, numeric(1), 1L)),
    bin_data_tuples
  )
}

# Corresponds to MassGrid.bin_track_mzs: aggregate m/z tuples of all samples into centroid bins.
MassGrid_bin_track_mzs <- function(self, tl, reference_id = NULL) {
  # The reference_id is not involved in the calculation for the time being, just like Python, only the interface is retained.
  invisible(reference_id)
  if (length(tl) == 0L) return(list())
  tolerance_ratio <- 1e-6 * self$experiment$parameters$mz_tolerance_ppm
  bins_of_bins <- list()
  current <- list(tl[[1L]])

  # First, a coarse bin is formed based on whether the adjacent m/z falls within the ppm tolerance.
  if (length(tl) > 1L) {
    for (ii in 2:length(tl)) {
      delta <- tl[[ii]][[1L]] - tl[[ii - 1L]][[1L]]
      if (delta < tolerance_ratio * tl[[ii - 1L]][[1L]]) {
        current[[length(current) + 1L]] <- tl[[ii]]
      } else {
        bins_of_bins[[length(bins_of_bins) + 1L]] <- current
        current <- list(tl[[ii]])
      }
    }
  }
  bins_of_bins[[length(bins_of_bins) + 1L]] <- current

  # When the range does not exceed double the tolerance, the median is taken directly, otherwise seed clustering is called to continue splitting.
  good_bins <- list()
  for (bin_data in bins_of_bins) {
    mz_values <- vapply(bin_data, `[[`, numeric(1), 1L)
    mz_range <- tail(mz_values, 1L) - mz_values[[1L]]
    mz_tolerance <- mz_values[[1L]] * tolerance_ratio
    clusters <- if (mz_range < mz_tolerance * 2) {
      list(bin_data)
    } else {
      nn_cluster_by_mz_seeds(bin_data, mz_tolerance)
    }
    for (cluster in clusters) {
      good_bins[[length(good_bins) + 1L]] <-
        MassGrid_bin_track_mzs__get_bin(cluster)
    }
  }
  good_bins
}

# Corresponds to MassGrid.join: Python’s original function is pass, so R remains side-effect-free.
MassGrid_join <- function(self, M2) {
  invisible(self)
  invisible(M2)
  NULL
}

# -----------------------------------------------------------------------------
# Python class CompositeMap: 20 methods
# -----------------------------------------------------------------------------

# Corresponds to CompositeMap.__init__: Create experimental composite map.
CompositeMap__init__ <- function(experiment) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("CompositeMap", "environment")
  self$experiment <- experiment
  self$`_number_of_samples_` <- experiment$number_of_samples

  # Sample column order strictly follows valid_sample_ids.
  self$list_sample_names <- unname(vapply(
    experiment$valid_sample_ids,
    function(sample_id) {
      .constructors_get(
        .constructors_registry_get(experiment$sample_registry, sample_id),
        "name"
      )
    },
    character(1)
  ))

  # reference sample determines the RT coordinate system and the maximum reference retention time.
  self$reference_sample_instance <-
    CompositeMap_get_reference_sample_instance(
      self, experiment$reference_sample_id
    )
  self$reference_sample <- self$reference_sample_instance
  self$rt_length <- experiment$number_scans
  self$dict_scan_rtime <- CompositeMap_get_reference_rtimes(
    self, self$rt_length
  )
  self$max_ref_rtime <- unname(
    self$dict_scan_rtime[[as.character(self$rt_length - 1L)]]
  )

  # Initialize the MassGrid, feature and composite track containers.
  self$MassGrid <- NULL
  self$FeatureTable <- NULL
  self$FeatureList <- list()
  self$`_mz_landmarks_` <- integer()
  self$good_reference_landmark_peaks <- list()
  self$composite_mass_tracks <- list()
  self
}

# Corresponds to CompositeMap.get_reference_sample_instance: wraps reference sample and loads mass tracks.
CompositeMap_get_reference_sample_instance <- function(self,
                                                       reference_sample_id) {
  registry <- .constructors_registry_get(
    self$experiment$sample_registry, reference_sample_id
  )
  sample <- .constructors_make_simple_sample(
    registry,
    self$experiment,
    self$experiment$database_mode,
    self$experiment$mode,
    is_reference = TRUE
  )
  sample$list_mass_tracks <- .constructors_get_mass_tracks(sample)
  sample
}

# Corresponds to CompositeMap.get_reference_rtimes: interpolates and extrapolates reference scan to retention time.
CompositeMap_get_reference_rtimes <- function(self, rt_length) {
  x <- self$reference_sample$rt_numbers
  y <- self$reference_sample$list_retention_time
  new_x <- seq_len(rt_length) - 1L
  new_y <- .constructors_interp_extrapolate(x, y, new_x)
  stats::setNames(new_y, new_x)
}

# Corresponds to CompositeMap.construct_mass_grid: Select sample-by-sample or centroiding path by number of samples.
CompositeMap_construct_mass_grid <- function(self) {
  message("Constructing MassGrid, ...")
  mass_grid <- MassGrid__init__(self, self$experiment)
  if (self$`_number_of_samples_` <=
      self$experiment$parameters$project_sample_number_small) {
    message("Building Grid Sample Wise...")
    MassGrid_build_grid_sample_wise(mass_grid)
  } else {
    message("Building Grid by Centroiding...")
    MassGrid_build_grid_by_centroiding(mass_grid)
  }
  self$MassGrid <- mass_grid$MassGrid
  self$`_mz_landmarks_` <- mass_grid$`_mz_landmarks_`
  invisible(self)
}

# Corresponds to CompositeMap.mock_retention_alignment: Creates an empty RT mapping for non-reference sample.
CompositeMap_mock_rentention_alignment <- function(self) {
  if (length(self$experiment$all_samples) > 1L) {
    for (sample in self$experiment$all_samples[-1L]) {
      sample$rt_cal_dict <- numeric()
      sample$reverse_rt_cal_dict <- numeric()
    }
  }
  invisible(self)
}

# Read the mapping target for the 0-based sample position in experiment.mapping.
.constructors_mapping_get <- function(mapping, sample_index) {
  if (is.null(mapping) || length(mapping) == 0L) return(NULL)
  key <- as.character(sample_index)
  if (!is.null(names(mapping)) && key %in% names(mapping)) {
    return(as.integer(mapping[[key]]))
  }
  NULL
}

# Corresponds to CompositeMap.perform_index_alignment: align and build composite tracks according to RI index samples.
CompositeMap_perform_index_alignment <- function(self) {
  mz_rows <- seq_len(nrow(self$MassGrid)) - 1L
  mz_values <- self$MassGrid$mz
  base_track <- integer(self$rt_length)
  composite <- lapply(mz_rows, function(index) base_track)

  # Samples that do not appear as ordinary study samples in the mapping are index standards.
  index_samples <- list()
  for (sample_index in seq_along(self$experiment$all_samples) - 1L) {
    if (is.null(.constructors_mapping_get(
      self$experiment$mapping, sample_index
    ))) {
      index_samples[[length(index_samples) + 1L]] <-
        self$experiment$all_samples[[sample_index + 1L]]
    }
  }
  if (length(index_samples) == 0L) {
    stop("Index alignment requires at least one index sample.", call. = FALSE)
  }

  self$good_reference_landmark_peaks <- CompositeMap_set_RT_reference(
    self, self$experiment$parameters$cal_min_peak_height
  )
  master_index_sample <- index_samples[[1L]]

  # Each index sample uses the same set of reference landmark peaks to build a scan map.
  for (index_sample in index_samples) {
    index_tracks <- .constructors_get_mass_tracks(index_sample)
    candidate_landmarks <- vapply(
      self$good_reference_landmark_peaks,
      function(peak) {
        .constructors_grid_get(
          self$MassGrid, index_sample$name, peak$ref_id_num
        )
      },
      numeric(1)
    )
    good_peaks <- list()
    selected_reference <- list()
    for (ii in seq_along(candidate_landmarks)) {
      track_id <- candidate_landmarks[[ii]]
      if (is.na(track_id)) next
      track <- index_tracks[[as.integer(track_id) + 1L]]
      unique_peak <- quick_detect_unique_elution_peak(
        track$intensity,
        min_peak_height = self$experiment$parameters$cal_min_peak_height,
        min_fwhm = 3,
        min_prominence_threshold_ratio = 0.2
      )
      if (!is.null(unique_peak)) {
        unique_peak$ref_id_num <- as.integer(track_id)
        good_peaks[[length(good_peaks) + 1L]] <- unique_peak
        selected_reference[[length(selected_reference) + 1L]] <-
          self$good_reference_landmark_peaks[[ii]]
      }
    }
    message("\tgood_landmarks: ", index_sample$name, " ", length(good_peaks))

    # Python manually performs clean points, LOWESS, and bidirectional linear extrapolation.
    sample_scans <- index_sample$list_scan_numbers
    sample_bound <- max(sample_scans)
    right_end <- 1.1 * sample_bound
    calibration_pairs <- Map(
      function(sample_peak, reference_peak) {
        c(sample_peak$apex, reference_peak$apex)
      },
      good_peaks,
      selected_reference
    )
    cleaned <- clean_rt_calibration_points(calibration_pairs)
    xx <- c(
      rep(-0.1 * sample_bound, 3L),
      vapply(cleaned, `[[`, numeric(1), 1L),
      rep(right_end, 3L)
    )
    yy <- c(
      rep(-0.1 * sample_bound, 3L),
      vapply(cleaned, `[[`, numeric(1), 2L),
      rep(right_end, 3L)
    )
    fraction <- max(0.2, min(
      0.6 - 0.004 * (length(cleaned) - 50), 0.6
    ))
    predicted <- hacked_lowess(
      yy, xx, fraction, 3, sample_scans
    )
    reverse <- .constructors_interp_extrapolate(
      predicted,
      sample_scans,
      master_index_sample$list_scan_numbers
    )

    predicted <- as.integer(round(predicted))
    forward_keep <- sample_scans != predicted & predicted >= 0L &
      predicted <= max(master_index_sample$list_scan_numbers)
    index_sample$rt_cal_dict <- stats::setNames(
      predicted[forward_keep], sample_scans[forward_keep]
    )
    reverse <- as.integer(round(reverse))
    reverse_keep <- master_index_sample$list_scan_numbers != reverse &
      reverse >= 0L & reverse <= sample_bound
    index_sample$reverse_rt_cal_dict <- stats::setNames(
      reverse[reverse_keep],
      master_index_sample$list_scan_numbers[reverse_keep]
    )
  }

  # Ordinary study samples are first mapped to the corresponding index samples through the RI model, and then superimposed on the composite.
  for (sample_index in seq_along(self$experiment$all_samples) - 1L) {
    sample <- self$experiment$all_samples[[sample_index + 1L]]
    tracks <- .constructors_get_mass_tracks(sample)
    message("Aligning: ", sample$name)
    index_position <- .constructors_mapping_get(
      self$experiment$mapping, sample_index
    )
    if (is.null(index_position)) next
    index_sample <- self$experiment$all_samples[[index_position + 1L]]
    ri_model <- self$experiment$RI_models[[index_position + 1L]]
    reverse_model <- self$experiment$reverse_RI_models[[index_position + 1L]]
    retention_indices <- ri_model(sample$list_retention_time)
    reference_scans <- reverse_model(retention_indices)

    forward <- reverse <- numeric()
    for (scan_index in seq_along(reference_scans) - 1L) {
      reference_scan <- as.integer(reference_scans[[scan_index + 1L]])
      fallback <- max(index_sample$list_scan_numbers)
      forward <- .constructors_dict_set(
        forward,
        scan_index,
        .constructors_dict_get(index_sample$rt_cal_dict, reference_scan, fallback)
      )
      reverse <- .constructors_dict_set(
        reverse,
        scan_index,
        .constructors_dict_get(
          index_sample$reverse_rt_cal_dict, reference_scan, fallback
        )
      )
    }
    sample$rt_cal_dict <- forward
    sample$reverse_rt_cal_dict <- reverse

    if (!self$experiment$parameters$drop_unaligned_samples ||
        isTRUE(sample$is_rt_aligned)) {
      for (row_id in mz_rows) {
        track_id <- .constructors_grid_get(
          self$MassGrid, sample$name, row_id
        )
        if (!is.na(track_id)) {
          remapped <- remap_intensity_track(
            tracks[[as.integer(track_id) + 1L]]$intensity,
            base_track,
            sample$rt_cal_dict
          )
          composite[[row_id + 1L]] <- composite[[row_id + 1L]] + remapped
        }
      }
    }
  }

  self$composite_mass_tracks <- lapply(seq_along(mz_rows), function(ii) {
    list(
      id_number = mz_rows[[ii]],
      mz = mz_values[[ii]],
      intensity = composite[[ii]]
    )
  })
  names(self$composite_mass_tracks) <- as.character(mz_rows)
  invisible(self)
}

# Corresponds to CompositeMap.build_composite_tracks_GC: GC path directly executes index alignment.
CompositeMap_build_composite_tracks_GC <- function(self) {
  CompositeMap_perform_index_alignment(self)
}

# Corresponds to CompositeMap.START inline def __similarity.
CompositeMap_START__similarity <- function(reference_peaks_per_sample,
                                           sample1,
                                           sample2) {
  indices1 <- unique(vapply(
    reference_peaks_per_sample[[sample1$name]], `[[`, numeric(1), "index"
  ))
  indices2 <- unique(vapply(
    reference_peaks_per_sample[[sample2$name]], `[[`, numeric(1), "index"
  ))
  union_indices <- union(indices1, indices2)

  # Python produces a divide-by-zero error on an empty union; R also explicitly rejects this undefined similarity.
  if (length(union_indices) == 0L) {
    stop("START similarity is undefined for two empty peak sets.", call. = FALSE)
  }
  length(intersect(indices1, indices2)) / length(union_indices)
}

# Corresponds to CompositeMap.START inline def __cost: minimum spanning tree uses 1-similarity.
CompositeMap_START__cost <- function(reference_peaks_per_sample,
                                     sample1,
                                     sample2) {
  1 - CompositeMap_START__similarity(
    reference_peaks_per_sample, sample1, sample2
  )
}

# Corresponds to CompositeMap.START embedded def __pairwise_cost.
CompositeMap_START__pairwise_cost <- function(samples,
                                              reference_peaks_per_sample) {
  vapply(samples, function(sample1) {
    vapply(samples, function(sample2) {
      CompositeMap_START__cost(
        reference_peaks_per_sample, sample1, sample2
      )
    }, numeric(1))
  }, numeric(length(samples)))
}

# Corresponds to CompositeMap.START embedded def __pairwise_similarity.
CompositeMap_START__pairwise_similarity <- function(samples,
                                                    reference_peaks_per_sample) {
  vapply(samples, function(sample1) {
    vapply(samples, function(sample2) {
      CompositeMap_START__similarity(
        reference_peaks_per_sample, sample1, sample2
      )
    }, numeric(1))
  }, numeric(length(samples)))
}

# Corresponds to CompositeMap.START embedded def __distance_to_graph.
CompositeMap_START__distance_to_graph <- function(distance_matrix) {
  graph <- as.matrix(distance_matrix)
  diag(graph) <- 0

  # networkx.from_numpy_array only interprets non-zero off-diagonal elements as edges.
  graph[graph == 0] <- NA_real_
  diag(graph) <- 0
  graph
}

# Replicate the core results of networkx.minimum_spanning_tree using Prim's algorithm.
.constructors_minimum_spanning_tree <- function(graph) {
  n <- nrow(graph)
  if (n <= 1L) return(graph)
  tree <- matrix(NA_real_, n, n)
  diag(tree) <- 0
  selected <- FALSE | seq_len(n) == 1L

  # Each round selects the minimum weight edge connecting the set to unconnected nodes.
  while (sum(selected) < n) {
    best_weight <- Inf
    best_from <- best_to <- NA_integer_
    for (from in which(selected)) {
      for (to in which(!selected)) {
        weight <- graph[from, to]
        if (!is.na(weight) && weight < best_weight) {
          best_weight <- weight
          best_from <- from
          best_to <- to
        }
      }
    }
    if (is.na(best_to)) stop("START distance graph is disconnected.", call. = FALSE)
    tree[best_from, best_to] <- tree[best_to, best_from] <- best_weight
    selected[[best_to]] <- TRUE
  }
  tree
}

# Computes the shortest distance in an unweighted graph, for use by the networkx.center equivalent implementation.
.constructors_graph_distances <- function(graph, start) {
  n <- nrow(graph)
  distance <- rep(Inf, n)
  distance[[start]] <- 0
  queue <- start
  while (length(queue) > 0L) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    neighbours <- which(!is.na(graph[node, ]) & seq_len(n) != node)
    for (neighbour in neighbours) {
      if (is.infinite(distance[[neighbour]])) {
        distance[[neighbour]] <- distance[[node]] + 1L
        queue <- c(queue, neighbour)
      }
    }
  }
  distance
}

# Corresponds to CompositeMap.START embedded def __find_graph_root.
CompositeMap_START__find_graph_root <- function(self,
                                                distance_graph,
                                                reference_peaks_per_sample) {
  eccentricity <- vapply(seq_len(nrow(distance_graph)), function(node) {
    max(.constructors_graph_distances(distance_graph, node))
  }, numeric(1))
  root <- which.min(eccentricity)

  # Python will align root with itself and write forward and reverse mappings.
  alignment <- CompositeMap_START__align_pair(
    self,
    self$experiment$all_samples[[root]],
    self$experiment$all_samples[[root]],
    reference_peaks_per_sample
  )
  self$experiment$all_samples[[root]]$rt_cal_dict <- alignment[[1L]]
  self$experiment$all_samples[[root]]$reverse_rt_cal_dict <- alignment[[2L]]
  root - 1L
}

# Corresponds to CompositeMap.START inline def __pairwise_traverse: Returns the only 0-based path on the tree.
CompositeMap_START__pairwise_traverse <- function(distance_graph,
                                                  root,
                                                  target) {
  root_r <- as.integer(root) + 1L
  target_r <- as.integer(target) + 1L
  n <- nrow(distance_graph)
  parent <- rep(NA_integer_, n)
  visited <- rep(FALSE, n)
  queue <- root_r
  visited[[root_r]] <- TRUE

  # BFS finds the unique shortest path in the tree.
  while (length(queue) > 0L && !visited[[target_r]]) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    neighbours <- which(!is.na(distance_graph[node, ]) & seq_len(n) != node)
    for (neighbour in neighbours) {
      if (!visited[[neighbour]]) {
        visited[[neighbour]] <- TRUE
        parent[[neighbour]] <- node
        queue <- c(queue, neighbour)
      }
    }
  }
  if (!visited[[target_r]]) return(NULL)

  path <- target_r
  while (path[[1L]] != root_r) path <- c(parent[[path[[1L]]]], path)
  as.integer(path - 1L)
}

# Corresponds to CompositeMap.START embedded def __align_pair.
CompositeMap_START__align_pair <- function(self,
                                           sample1,
                                           sample2,
                                           reference_peaks_per_sample) {
  message("Aligning: ", sample1$name, " to ", sample2$name)
  peaks1 <- reference_peaks_per_sample[[sample1$name]]
  peaks2 <- reference_peaks_per_sample[[sample2$name]]
  mzs1 <- unique(vapply(peaks1, `[[`, numeric(1), "mz"))
  mzs2 <- unique(vapply(peaks2, `[[`, numeric(1), "mz"))
  shared <- intersect(mzs1, mzs2)
  message(
    "\tPeaks - Shared / Sample 1 / Sample 2: ",
    length(shared), " / ", length(mzs1), " / ", length(mzs2)
  )

  # Record the apex of two samples for each shared three-digit decimal m/z.
  reference_pairs <- list()
  for (mz in shared) {
    apex1 <- peaks1[[which(vapply(peaks1, `[[`, numeric(1), "mz") == mz)[[1L]]]]$apex
    apex2 <- peaks2[[which(vapply(peaks2, `[[`, numeric(1), "mz") == mz)[[1L]]]]$apex
    reference_pairs[[as.character(mz)]] <- c(apex1, apex2)
  }
  cleaned <- clean_rt_calibration_points(unname(reference_pairs))
  reference_scans <- sample1$rt_numbers
  sample_scans <- sample2$rt_numbers
  reference_bound <- max(reference_scans)
  sample_bound <- max(sample_scans)
  right_end <- 1.1 * sample_bound
  xx <- c(
    rep(-0.1 * sample_bound, 3L),
    vapply(cleaned, `[[`, numeric(1), 1L),
    rep(right_end, 3L)
  )
  yy <- c(
    rep(-0.1 * sample_bound, 3L),
    vapply(cleaned, `[[`, numeric(1), 2L),
    rep(right_end, 3L)
  )
  predicted <- hacked_lowess(yy, xx, 0.5, 3, sample_scans)
  reverse <- .constructors_interp_extrapolate(
    predicted, sample_scans, reference_scans
  )

  # START The original implementation retains LOWESS floating point results and does not round here.
  keep_forward <- sample_scans != predicted & predicted >= 0 &
    predicted <= reference_bound
  forward_dict <- stats::setNames(
    predicted[keep_forward], sample_scans[keep_forward]
  )
  keep_reverse <- reference_scans != reverse & reverse >= 0 &
    reverse <= sample_bound
  reverse_dict <- stats::setNames(
    reverse[keep_reverse], reference_scans[keep_reverse]
  )
  list(forward_dict, reverse_dict)
}

# Corresponds to CompositeMap.START embedded def __align.
CompositeMap_START__align <- function(self,
                                      path,
                                      reference_peaks_per_sample) {
  calibrated_domain <- self$experiment$all_samples[[path[[1L]] + 1L]]$rt_numbers
  target_sample <- self$experiment$all_samples[[path[[1L]] + 1L]]

  # The RT maps for each pair of samples are combined sequentially along the spanning tree path.
  if (length(path) > 1L) {
    for (ii in seq_len(length(path) - 1L)) {
      sample1 <- self$experiment$all_samples[[path[[ii]] + 1L]]
      sample2 <- self$experiment$all_samples[[path[[ii + 1L]] + 1L]]
      alignment <- CompositeMap_START__align_pair(
        self, sample1, sample2, reference_peaks_per_sample
      )
      calibrated_domain <- vapply(calibrated_domain, function(scan) {
        .constructors_dict_get(alignment[[1L]], scan, scan)
      }, numeric(1))
      target_sample <- sample2
    }
  }
  calibrated_domain <- as.integer(round(calibrated_domain))
  source_scans <- self$experiment$all_samples[[path[[1L]] + 1L]]$rt_numbers
  forward <- stats::setNames(calibrated_domain, source_scans)
  reverse <- stats::setNames(source_scans, calibrated_domain)
  target_sample$rt_cal_dict <- forward
  target_sample$reverse_rt_cal_dict <- reverse
  target_sample$is_rt_aligned <- TRUE
  list(forward, reverse)
}

# Corresponds to CompositeMap.START: Experimental spanning-tree RT alignment full path.
CompositeMap_START <- function(self) {
  cal_min_peak_height <- self$experiment$parameters$cal_min_peak_height
  min_selectivity <- 0.99
  selectivities <- calculate_selectivity(
    self$MassGrid$mz,
    self$experiment$parameters$mz_tolerance_ppm
  )

  # Create an inverted index of (sample_name, track_id) into a MassGrid 0-based row number.
  mass_grid_lookup <- list()
  for (row in seq_len(nrow(self$MassGrid))) {
    for (sample_name in names(self$MassGrid)[-1L]) {
      track_id <- self$MassGrid[[sample_name]][[row]]
      if (!is.na(track_id)) {
        key <- paste(sample_name, track_id, sep = "\r")
        mass_grid_lookup[[key]] <- row - 1L
      }
    }
  }

  # Only landmarks with high m/z selectivity and unique prominent peaks are retained for each sample.
  reference_peaks_per_sample <- list()
  for (sample in self$experiment$all_samples) {
    sample_peaks <- list()
    tracks <- .constructors_get_mass_tracks(sample)
    for (track_index in seq_along(tracks) - 1L) {
      key <- paste(sample$name, track_index, sep = "\r")
      mapped_index <- mass_grid_lookup[[key]]
      if (is.null(mapped_index) || selectivities[[mapped_index + 1L]] <=
          min_selectivity) next
      peak <- quick_detect_unique_elution_peak(
        tracks[[track_index + 1L]]$intensity,
        min_peak_height = cal_min_peak_height,
        min_fwhm = 3,
        min_prominence_threshold_ratio = 0.2
      )
      if (!is.null(peak)) {
        peak$mz <- round(tracks[[track_index + 1L]]$mz, 3)
        peak$index <- mapped_index
        sample_peaks[[length(sample_peaks) + 1L]] <- peak
      }
    }
    reference_peaks_per_sample[[sample$name]] <- sample_peaks
  }

  # After the similarity is converted to distance, a minimum spanning tree is established, and the calibration is propagated from the center of the tree to all nodes.
  distances <- CompositeMap_START__pairwise_cost(
    self$experiment$all_samples,
    reference_peaks_per_sample
  )
  graph <- CompositeMap_START__distance_to_graph(distances)
  tree <- .constructors_minimum_spanning_tree(graph)
  root <- CompositeMap_START__find_graph_root(
    self, tree, reference_peaks_per_sample
  )
  for (node in seq_len(nrow(tree)) - 1L) {
    if (node == root) next
    path <- CompositeMap_START__pairwise_traverse(tree, root, node)
    CompositeMap_START__align(self, path, reference_peaks_per_sample)
  }

  # Overlay all mass tracks using the map obtained for each sample.
  mz_rows <- seq_len(nrow(self$MassGrid)) - 1L
  base_track <- integer(self$rt_length)
  composite <- lapply(mz_rows, function(index) base_track)
  for (sample in self$experiment$all_samples) {
    sample$is_rt_aligned <- TRUE
    tracks <- .constructors_get_mass_tracks(sample)
    for (row_id in mz_rows) {
      track_id <- .constructors_grid_get(self$MassGrid, sample$name, row_id)
      if (!is.na(track_id)) {
        composite[[row_id + 1L]] <- composite[[row_id + 1L]] +
          remap_intensity_track(
            tracks[[as.integer(track_id) + 1L]]$intensity,
            base_track,
            sample$rt_cal_dict
          )
      }
    }
  }
  self$composite_mass_tracks <- lapply(seq_along(mz_rows), function(ii) {
    list(
      id_number = mz_rows[[ii]],
      mz = self$MassGrid$mz[[ii]],
      intensity = composite[[ii]]
    )
  })
  names(self$composite_mass_tracks) <- as.character(mz_rows)
  invisible(self)
}

# Corresponds to CompositeMap.build_composite_tracks: RT calibration and track overlay for regular LC paths.
CompositeMap_build_composite_tracks <- function(self) {
  message("\nBuilding composite mass tracks and calibrating retention time ...\n")
  parameters <- self$experiment$parameters
  cal_min_peak_height <- parameters$cal_min_peak_height
  min_peak_number <- parameters$peak_number_rt_calibration
  number_iterations <- parameters$num_lowess_iterations
  max_shift <- if (is.null(parameters$max_retention_shift)) {
    Inf
  } else {
    parameters$max_retention_shift
  }
  self$good_reference_landmark_peaks <- CompositeMap_set_RT_reference(
    self, cal_min_peak_height
  )

  # Each MassGrid row establishes a composite base track with int64 semantics.
  mz_rows <- seq_len(nrow(self$MassGrid)) - 1L
  base_track <- integer(self$rt_length)
  composite <- lapply(mz_rows, function(index) base_track)

  # The debug mode first exports the reference landmark table like Python.
  if (isTRUE(parameters$debug_rtime_align)) {
    CompositeMap_export_reference_sample(self)
  }

  # Calibrate RT sample by sample, and then add remapped intensity to the corresponding composite row.
  for (sample in self$experiment$all_samples) {
    message("   ", sample$name)
    tracks <- .constructors_get_mass_tracks(sample)
    if (isTRUE(sample$is_reference)) {
      message(
        "\t\tgood_reference_landmark_peaks: ",
        length(self$good_reference_landmark_peaks)
      )
    } else if (isTRUE(parameters$rt_align_on)) {
      calibration_function <- if (isTRUE(parameters$debug_rtime_align)) {
        rt_lowess_calibration_debug
      } else {
        rt_lowess_calibration
      }
      CompositeMap_calibrate_sample_RT(
        self,
        sample,
        tracks,
        calibration_fuction = calibration_function,
        cal_min_peak_height = cal_min_peak_height,
        MIN_PEAK_NUM = min_peak_number,
        MAX_RETENTION_SHIFT = max_shift,
        NUM_ITERATIONS = number_iterations
      )
    }

    if (!isTRUE(parameters$drop_unaligned_samples) ||
        isTRUE(sample$is_rt_aligned)) {
      for (row_id in mz_rows) {
        track_id <- .constructors_grid_get(
          self$MassGrid, sample$name, row_id
        )
        if (!is.na(track_id)) {
          remapped <- remap_intensity_track(
            tracks[[as.integer(track_id) + 1L]]$intensity,
            base_track,
            sample$rt_cal_dict
          )
          composite[[row_id + 1L]] <- composite[[row_id + 1L]] + remapped
        }
      }
    }
  }

  self$composite_mass_tracks <- lapply(seq_along(mz_rows), function(ii) {
    list(
      id_number = mz_rows[[ii]],
      mz = self$MassGrid$mz[[ii]],
      intensity = composite[[ii]]
    )
  })
  names(self$composite_mass_tracks) <- as.character(mz_rows)
  invisible(self)
}

# Corresponds to CompositeMap.calibrate_sample_RT_by_standards: Python original version is pass.
CompositeMap_calibrate_sample_RT_by_standards <- function(self, sample) {
  invisible(self)
  invisible(sample)
  NULL
}

# Corresponds to CompositeMap.calibrate_sample_RT: Calibrate a sample with unique landmark peaks.
CompositeMap_calibrate_sample_RT <- function(
    self,
    sample,
    list_mass_tracks,
    calibration_fuction = rt_lowess_calibration,
    cal_min_peak_height = 100000,
    MIN_PEAK_NUM = 15,
    MAX_RETENTION_SHIFT = Inf,
    NUM_ITERATIONS = 3) {
  # Press the MassGrid row of the reference landmark to find the corresponding track ID of the sample.
  candidate_landmarks <- vapply(
    self$good_reference_landmark_peaks,
    function(reference_peak) {
      .constructors_grid_get(
        self$MassGrid, sample$name, reference_peak$ref_id_num
      )
    },
    numeric(1)
  )
  good_peaks <- list()
  selected_reference <- list()

  # Only landmarks with unique prominent peaks and RT shifts that do not exceed the threshold are retained.
  for (ii in seq_along(candidate_landmarks)) {
    track_id <- candidate_landmarks[[ii]]
    if (is.na(track_id)) next
    track <- list_mass_tracks[[as.integer(track_id) + 1L]]
    unique_peak <- quick_detect_unique_elution_peak(
      track$intensity,
      min_peak_height = cal_min_peak_height,
      min_fwhm = 3,
      min_prominence_threshold_ratio = 0.2
    )
    if (!is.null(unique_peak)) {
      scan_delta <- unique_peak$apex -
        self$good_reference_landmark_peaks[[ii]]$apex
      if (abs(scan_delta) < MAX_RETENTION_SHIFT) {
        unique_peak$ref_id_num <- as.integer(track_id)
        good_peaks[[length(good_peaks) + 1L]] <- unique_peak
        selected_reference[[length(selected_reference) + 1L]] <-
          self$good_reference_landmark_peaks[[ii]]
      }
    }
  }
  message("\tgood_landmark_peaks: ", length(good_peaks))
  sample$rt_landmarks <- if (length(good_peaks) == 0L) integer() else {
    vapply(good_peaks, `[[`, numeric(1), "apex")
  }

  # Python uses strict > MIN_PEAK_NUM; the threshold is met but not exceeded and still not calibrated.
  if (length(good_peaks) > MIN_PEAK_NUM) {
    calibration <- tryCatch(
      calibration_fuction(
        good_peaks,
        selected_reference,
        sample$rt_numbers,
        self$reference_sample$rt_numbers,
        NUM_ITERATIONS,
        sample$name,
        self$experiment$parameters$outdir
      ),
      error = function(error) NULL
    )
    if (!is.null(calibration)) {
      sample$rt_cal_dict <- calibration[[1L]]
      sample$reverse_rt_cal_dict <- calibration[[2L]]
      sample$is_rt_aligned <- TRUE
    }
  }

  # When calibration fails, an empty map is retained so that the sample can still be stacked at the original scan position.
  if (!isTRUE(sample$is_rt_aligned)) {
    sample$rt_cal_dict <- numeric()
    sample$reverse_rt_cal_dict <- numeric()
    message(
      "    ~warning~ Failure in retention time alignment (",
      length(good_peaks), "); ", sample$name, "."
    )
  }
  invisible(sample)
}

# Corresponds to CompositeMap.set_RT_reference: select the only peak with high selectivity in reference sample.
CompositeMap_set_RT_reference <- function(self,
                                          cal_peak_intensity_threshold = 100000) {
  landmarks <- as.integer(self$`_mz_landmarks_`)
  if (length(landmarks) == 0L) {
    self$reference_sample$rt_landmarks <- integer()
    return(list())
  }
  selectivities <- calculate_selectivity(
    self$MassGrid$mz[landmarks],
    self$experiment$parameters$mz_tolerance_ppm
  )
  reference_tracks <- self$reference_sample$list_mass_tracks
  good_reference <- list()

  # Python `if ref_ii` excludes track ID 0; this original behavior is explicitly retained here.
  for (ii in seq_along(landmarks)) {
    if (selectivities[[ii]] <= 0.99) next
    row_id <- landmarks[[ii]] - 1L
    track_id <- .constructors_grid_get(
      self$MassGrid, self$reference_sample$name, row_id
    )
    if (is.na(track_id) || track_id == 0) next
    track <- reference_tracks[[as.integer(track_id) + 1L]]
    unique_peak <- quick_detect_unique_elution_peak(
      track$intensity,
      min_peak_height = cal_peak_intensity_threshold,
      min_fwhm = 3,
      min_prominence_threshold_ratio = 0.2
    )
    if (!is.null(unique_peak)) {
      unique_peak$ref_id_num <- row_id
      good_reference[[length(good_reference) + 1L]] <- unique_peak
    }
  }
  self$reference_sample$rt_landmarks <- if (length(good_reference) == 0L) {
    integer()
  } else {
    vapply(good_reference, `[[`, numeric(1), "apex")
  }
  good_reference
}

# Convert named peak lists to a data.frame, preserving the order of first occurrence of all fields.
.constructors_records_to_data_frame <- function(records) {
  if (length(records) == 0L) return(data.frame())
  columns <- unique(unlist(lapply(records, names), use.names = FALSE))
  result <- lapply(columns, function(column) {
    values <- lapply(records, function(record) record[[column]])
    # The feature fields of constructors are all length-scalars; missing fields are NA.
    vapply(values, function(value) {
      if (is.null(value) || length(value) == 0L) NA else value[[1L]]
    }, FUN.VALUE = if (all(vapply(values, function(value) {
      is.null(value) || is.character(value)
    }, logical(1)))) character(1) else numeric(1))
  })
  names(result) <- columns
  as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
}

# Corresponds to CompositeMap.global_peak_detection: detects experimental features on composite tracks.
CompositeMap_global_peak_detection <- function(self) {
  message(
    "\nPeak detection on ", length(self$composite_mass_tracks),
    " composite mass tracks, ...\n"
  )
  self$FeatureList <- batch_deep_detect_elution_peaks(
    unname(self$composite_mass_tracks),
    self$experiment$number_scans,
    self$experiment$parameters
  )

  # Assign F0, F1... to each peak and convert scan boundaries to reference rtime.
  for (ii in seq_along(self$FeatureList)) {
    peak <- self$FeatureList[[ii]]
    peak$id_number <- paste0("F", ii - 1L)
    apex_time <- .constructors_dict_get(
      self$dict_scan_rtime, peak$apex, NULL
    )
    if (is.null(apex_time)) {
      apex_time <- self$max_ref_rtime
      message("Feature rtime out of bound - ", peak$id_number, " ", peak$apex)
    }
    peak$rtime <- apex_time
    left_time <- .constructors_dict_get(
      self$dict_scan_rtime, peak$left_base, NULL
    )
    right_time <- .constructors_dict_get(
      self$dict_scan_rtime, peak$right_base, NULL
    )
    if (!is.null(left_time) && !is.null(right_time)) {
      peak$rtime_left_base <- left_time
      peak$rtime_right_base <- right_time
    } else {
      message("Feature rtime out of bound on ", peak$id_number)
    }
    self$FeatureList[[ii]] <- peak
  }
  CompositeMap_generate_feature_table(self)
  invisible(self)
}

# Corresponds to CompositeMap.get_peak_area_sum: summation of intensity within a closed interval.
CompositeMap_get_peak_area_sum <- function(self,
                                           track_intensity,
                                           left_base,
                                           right_base) {
  invisible(self)
  left_base <- as.integer(left_base)
  right_base <- as.integer(right_base)
  sum(.constructors_python_slice(
    track_intensity, left_base, right_base + 1L
  ))
}

# Corresponds to CompositeMap.get_peak_area_auc: size=2 maximum filter, then sum and int().
CompositeMap_get_peak_area_auc <- function(self,
                                           track_intensity,
                                           left_base,
                                           right_base) {
  invisible(self)
  interval <- .constructors_python_slice(
    track_intensity,
    as.integer(left_base),
    as.integer(right_base) + 1L
  )
  as.integer(trunc(sum(.constructors_maximum_filter_size2(interval))))
}

# Corresponds to CompositeMap.get_peak_area_gaussian: Gaussian fully integrated and truncated by Python int().
CompositeMap_get_peak_area_gaussian <- function(self,
                                                track_intensity,
                                                left_base,
                                                right_base) {
  invisible(self)
  area <- get_gaussian_peakarea_on_intensity_list(
    track_intensity, left_base, right_base
  )
  if (exists(".python_int", mode = "function")) return(.python_int(area))
  trunc(area)
}

# Corresponds to CompositeMap.get_DIMS_feature_table: each mass track uses the maximum intensity as the feature.
CompositeMap_get_DIMS_feature_table <- function(self) {
  self$FeatureList <- lapply(seq_len(nrow(self$MassGrid)), function(row) {
    list(
      id_number = paste0("F", row - 1L),
      parent_masstrack_id = row - 1L,
      mz = self$MassGrid$mz[[row]],
      intensity = 0,
      apex = 0,
      left_base = 0,
      right_base = 0,
      rtime = 0,
      rtime_left_base = 0,
      rtime_right_base = 0,
      peak_area = 0,
      cSelectivity = 0,
      goodness_fitting = 0,
      snr = 0
    )
  })
  feature_table <- .constructors_records_to_data_frame(self$FeatureList)

  # DIMS does not use peak boundaries, and each sample column directly takes the maximum intensity of the corresponding mass track.
  for (sample in self$experiment$all_samples) {
    tracks <- .constructors_get_mass_tracks(sample)
    feature_table[[sample$name]] <- vapply(self$FeatureList, function(feature) {
      track_id <- .constructors_grid_get(
        self$MassGrid, sample$name, feature$parent_masstrack_id
      )
      if (is.na(track_id)) 0 else {
        max(tracks[[as.integer(track_id) + 1L]]$intensity)
      }
    }, numeric(1))
  }
  self$FeatureTable <- feature_table
  invisible(self)
}

# Corresponds to CompositeMap.generate_feature_table: Select the peak area algorithm according to parameters and fill the sample column.
CompositeMap_generate_feature_table <- function(self) {
  area_methods <- list(
    auc = CompositeMap_get_peak_area_auc,
    sum = CompositeMap_get_peak_area_sum,
    gauss = CompositeMap_get_peak_area_gaussian
  )
  area_function <- area_methods[[self$experiment$parameters$peak_area]]
  if (is.null(area_function)) stop("Unknown peak_area method.", call. = FALSE)
  feature_table <- .constructors_records_to_data_frame(self$FeatureList)

  # Python will skip samples that are not aligned even when drop_unaligned_samples is enabled.
  for (sample in self$experiment$all_samples) {
    if (!isTRUE(self$experiment$parameters$drop_unaligned_samples) ||
        isTRUE(sample$is_rt_aligned)) {
      feature_table[[sample$name]] <- CompositeMap_extract_features_per_sample(
        self, sample, area_function
      )
    }
  }
  self$FeatureTable <- feature_table
  invisible(self)
}

# Corresponds to CompositeMap.extract_features_per_sample: extract the sample area according to the experiment feature boundary.
CompositeMap_extract_features_per_sample <- function(self,
                                                     sample,
                                                     peak_area_function) {
  tracks <- .constructors_get_mass_tracks(sample)
  vapply(self$FeatureList, function(feature) {
    track_id <- .constructors_grid_get(
      self$MassGrid, sample$name, feature$parent_masstrack_id
    )
    if (is.na(track_id)) return(0)
    track <- tracks[[as.integer(track_id) + 1L]]
    left_base <- .constructors_dict_get(
      sample$reverse_rt_cal_dict, feature$left_base, feature$left_base
    )
    right_base <- .constructors_dict_get(
      sample$reverse_rt_cal_dict, feature$right_base, feature$right_base
    )
    peak_area_function(self, track$intensity, left_base, right_base)
  }, numeric(1))
}

# Corresponds to CompositeMap.export_reference_sample: write the mz/rtime CSV of the reference landmark.
CompositeMap_export_reference_sample <- function(self) {
  mz_landmarks <- vapply(
    self$good_reference_landmark_peaks,
    function(peak) self$MassGrid$mz[[peak$ref_id_num + 1L]],
    numeric(1)
  )
  rtime_landmarks <- vapply(
    self$good_reference_landmark_peaks,
    function(peak) {
      .constructors_dict_get(self$dict_scan_rtime, peak$apex, NA_real_)
    },
    numeric(1)
  )
  output_path <- file.path(
    self$experiment$parameters$outdir,
    "export",
    paste0(self$reference_sample$name, "_mz_rtime_landmarks.csv")
  )

  # Like Python csv.writer, write the mz and rtime column headers and data rows.
  utils::write.csv(
    data.frame(mz = mz_landmarks, rtime = rtime_landmarks),
    output_path,
    row.names = FALSE,
    quote = FALSE
  )
  invisible(output_path)
}

# Provide R constructor aliases with the same Python class names; no new algorithm functions are added.
MassGrid <- MassGrid__init__
CompositeMap <- CompositeMap__init__
