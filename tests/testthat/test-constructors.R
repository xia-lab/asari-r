# Each of the 37 defs in constructors.py must have an explicit R corresponding function.
test_that("all 37 Python constructors defs have explicit R counterparts", {
  expected_functions <- c(
    "MassGrid__init__",
    "MassGrid_build_grid_sample_wise",
    "MassGrid_build_grid_by_centroiding",
    "MassGrid__initiate_mass_grid",
    "MassGrid_add_sample",
    "MassGrid_bin_track_mzs",
    "MassGrid_join",
    "MassGrid_bin_track_mzs__get_bin",
    "CompositeMap__init__",
    "CompositeMap_get_reference_sample_instance",
    "CompositeMap_get_reference_rtimes",
    "CompositeMap_construct_mass_grid",
    "CompositeMap_mock_rentention_alignment",
    "CompositeMap_perform_index_alignment",
    "CompositeMap_build_composite_tracks_GC",
    "CompositeMap_START",
    "CompositeMap_build_composite_tracks",
    "CompositeMap_calibrate_sample_RT_by_standards",
    "CompositeMap_calibrate_sample_RT",
    "CompositeMap_set_RT_reference",
    "CompositeMap_global_peak_detection",
    "CompositeMap_get_peak_area_sum",
    "CompositeMap_get_peak_area_auc",
    "CompositeMap_get_peak_area_gaussian",
    "CompositeMap_get_DIMS_feature_table",
    "CompositeMap_generate_feature_table",
    "CompositeMap_extract_features_per_sample",
    "CompositeMap_export_reference_sample",
    "CompositeMap_START__similarity",
    "CompositeMap_START__cost",
    "CompositeMap_START__pairwise_cost",
    "CompositeMap_START__pairwise_similarity",
    "CompositeMap_START__distance_to_graph",
    "CompositeMap_START__find_graph_root",
    "CompositeMap_START__pairwise_traverse",
    "CompositeMap_START__align_pair",
    "CompositeMap_START__align"
  )

  # The number must be exactly 37, and each name resolves to a function.
  expect_length(expected_functions, 37L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# Creates a minimal mutable SimpleSample used by constructors tests.
make_constructor_sample <- function(name,
                                    mzs,
                                    intensities,
                                    is_reference = FALSE) {
  sample <- new.env(parent = emptyenv())
  class(sample) <- c("SimpleSample", "environment")
  sample$name <- name
  sample$is_reference <- is_reference
  sample$is_rt_aligned <- is_reference
  sample$rt_numbers <- seq_along(intensities[[1L]]) - 1L
  sample$list_scan_numbers <- sample$rt_numbers
  sample$list_retention_time <- sample$rt_numbers * 10
  sample$rt_cal_dict <- numeric()
  sample$reverse_rt_cal_dict <- numeric()
  sample$`_mz_landmarks_` <- integer()
  sample$list_mass_tracks <- lapply(seq_along(mzs), function(ii) {
    list(
      id_number = ii - 1L,
      mz = mzs[[ii]],
      intensity = intensities[[ii]]
    )
  })
  sample$track_mzs <- lapply(seq_along(mzs), function(ii) {
    c(mzs[[ii]], ii - 1L)
  })
  sample$get_masstracks_and_anchors <- function() sample$list_mass_tracks
  sample
}

# Create a minimal experimental environment that can be used by CompositeMap and MassGrid.
make_constructor_experiment <- function(samples) {
  experiment <- new.env(parent = emptyenv())
  sample_ids <- names(samples)
  experiment$sample_registry <- samples
  experiment$valid_sample_ids <- sample_ids
  experiment$reference_sample_id <- sample_ids[[1L]]
  experiment$number_of_samples <- length(samples)
  experiment$number_scans <- length(samples[[1L]]$rt_numbers)
  experiment$database_mode <- "memory"
  experiment$mode <- "pos"
  experiment$all_samples <- list()
  experiment$parameters <- list(
    mz_tolerance_ppm = 5,
    correction_tolerance_ppm = 1,
    project_sample_number_small = 10,
    cal_min_peak_height = 100000,
    peak_number_rt_calibration = 15,
    num_lowess_iterations = 3,
    max_retention_shift = NULL,
    debug_rtime_align = FALSE,
    rt_align_on = FALSE,
    drop_unaligned_samples = FALSE,
    peak_area = "sum",
    outdir = tempdir(),
    min_peak_height = 100000,
    min_peak_ratio = 0.001,
    min_timepoints = 6,
    min_intensity_threshold = 1000,
    wlen = 25,
    signal_noise_ratio = 2,
    gaussian_shape = 0.5,
    min_prominence_threshold = 33000
  )
  experiment
}

# Verify that the R environment constructor preserves Python object fields and mutable state.
test_that("MassGrid and CompositeMap constructors initialize Python fields", {
  reference <- make_constructor_sample(
    "reference", c(100, 101), list(1:5, 5:1), TRUE
  )
  experiment <- make_constructor_experiment(list(reference = reference))
  composite <- CompositeMap__init__(experiment)
  mass_grid <- MassGrid__init__(composite, experiment)

  expect_s3_class(composite, "CompositeMap")
  expect_s3_class(mass_grid, "MassGrid")
  expect_identical(composite$reference_sample, reference)
  expect_equal(composite$dict_scan_rtime, stats::setNames(0:4 * 10, 0:4))
  expect_equal(composite$max_ref_rtime, 40)
  expect_equal(mass_grid$list_sample_names, "reference")
})

# Python reference results: the first three m/z belong to a coarse bin, and 101 is binned separately.
test_that("MassGrid bin_track_mzs reproduces ppm centroid bins", {
  self <- new.env(parent = emptyenv())
  self$experiment <- new.env(parent = emptyenv())
  self$experiment$parameters <- list(mz_tolerance_ppm = 5)
  tuples <- list(
    c(100.0000, 0, 0),
    c(100.0002, 0, 1),
    c(100.0003, 1, 0),
    c(101.0000, 1, 1)
  )

  bins <- MassGrid_bin_track_mzs(self, tuples)

  expect_length(bins, 2L)
  expect_equal(bins[[1L]][[1L]], 100.0002)
  expect_length(bins[[1L]][[2L]], 3L)
  expect_equal(bins[[2L]][[1L]], 101)
})

# Verify that large sample centroid paths generate mz columns and 0-based track ID sample columns.
test_that("centroid grid aligns tracks from multiple samples", {
  sample1 <- make_constructor_sample(
    "s1", c(100, 101), list(1:5, 5:1), TRUE
  )
  sample2 <- make_constructor_sample(
    "s2", c(100.0001, 102), list(2:6, 6:2), FALSE
  )
  experiment <- make_constructor_experiment(list(s1 = sample1, s2 = sample2))
  composite <- new.env(parent = emptyenv())
  composite$reference_sample_instance <- sample1
  composite$max_ref_rtime <- 40
  composite$list_sample_names <- c("s1", "s2")
  composite$`_number_of_samples_` <- 2L
  mass_grid <- MassGrid__init__(composite, experiment)

  MassGrid_build_grid_by_centroiding(mass_grid)

  expect_equal(mass_grid$MassGrid$mz, c(100.00005, 101, 102))
  expect_equal(mass_grid$MassGrid$s1, c(0, 1, NA))
  expect_equal(mass_grid$MassGrid$s2, c(0, NA, 1))
  expect_length(experiment$all_samples, 2L)
})

# Verifying the small sample path will initialize the reference grid and add new samples through landmark mapping.
test_that("sample-wise grid initializes reference and adds every sample", {
  sample1 <- make_constructor_sample(
    "s1", c(100, 101), list(1:5, 5:1), TRUE
  )
  sample2 <- make_constructor_sample(
    "s2", c(100.0001, 102), list(2:6, 6:2), FALSE
  )
  experiment <- make_constructor_experiment(list(s1 = sample1, s2 = sample2))
  composite <- CompositeMap(experiment)

  CompositeMap_construct_mass_grid(composite)

  expect_equal(composite$MassGrid$mz, c(100.00005, 101, 102))
  expect_equal(composite$MassGrid$s1, c(0, 1, NA))
  expect_equal(composite$MassGrid$s2, c(0, NA, 1))
  expect_length(experiment$all_samples, 2L)
  expect_identical(MassGrid, MassGrid__init__)
  expect_identical(CompositeMap, CompositeMap__init__)
})

# Validation of SciPy interp1d for interval interpolation and extrapolation of both ends.
test_that("reference retention times interpolate and extrapolate", {
  self <- new.env(parent = emptyenv())
  self$reference_sample <- make_constructor_sample(
    "reference", 100, list(1:3), TRUE
  )
  self$reference_sample$rt_numbers <- c(0, 2, 4)
  self$reference_sample$list_retention_time <- c(0, 20, 40)

  result <- CompositeMap_get_reference_rtimes(self, 6)

  expect_equal(unname(result), c(0, 10, 20, 30, 40, 50))
  expect_equal(names(result), as.character(0:5))
})

# Python scipy.maximum_filter1d(size=2, constant) returns 1,5,5,4 for this vector.
test_that("peak area methods reproduce Python boundaries and maximum filter", {
  self <- new.env(parent = emptyenv())
  intensity <- c(1, 5, 2, 4)

  expect_equal(CompositeMap_get_peak_area_sum(self, intensity, 0, 3), 12)
  expect_equal(CompositeMap_get_peak_area_sum(self, intensity, 1, 2), 7)
  expect_equal(CompositeMap_get_peak_area_auc(self, intensity, 0, 3), 15L)
  expect_null(CompositeMap_calibrate_sample_RT_by_standards(self, NULL))
  expect_null(MassGrid_join(self, NULL))
})

# Verify set similarity, distance matrix, and tree path behavior of START nested defs.
test_that("START helper defs reproduce similarity and graph traversal", {
  sample1 <- new.env(parent = emptyenv())
  sample2 <- new.env(parent = emptyenv())
  sample1$name <- "s1"
  sample2$name <- "s2"
  peaks <- list(
    s1 = list(list(index = 1), list(index = 2), list(index = 3)),
    s2 = list(list(index = 2), list(index = 3), list(index = 4))
  )

  expect_equal(
    CompositeMap_START__similarity(peaks, sample1, sample2),
    0.5
  )
  expect_equal(CompositeMap_START__cost(peaks, sample1, sample2), 0.5)
  expect_equal(
    CompositeMap_START__pairwise_similarity(list(sample1, sample2), peaks),
    matrix(c(1, 0.5, 0.5, 1), 2)
  )

  tree <- matrix(NA_real_, 4, 4)
  diag(tree) <- 0
  tree[1, 2] <- tree[2, 1] <- 1
  tree[2, 3] <- tree[3, 2] <- 1
  tree[3, 4] <- tree[4, 3] <- 1
  expect_equal(
    CompositeMap_START__pairwise_traverse(tree, 0, 3),
    0:3
  )
})

# Verify that the regular composite path overlays sample intensity according to the MassGrid correspondence.
test_that("build_composite_tracks sums aligned sample mass tracks", {
  reference <- make_constructor_sample(
    "reference", 100, list(as.integer(c(0, 10, 20, 10, 0))), TRUE
  )
  sample <- make_constructor_sample(
    "sample", 100.0001, list(as.integer(c(0, 5, 15, 5, 0))), FALSE
  )
  sample$is_rt_aligned <- FALSE
  experiment <- make_constructor_experiment(
    list(reference = reference, sample = sample)
  )
  experiment$all_samples <- list(reference, sample)
  composite <- CompositeMap__init__(experiment)
  composite$MassGrid <- data.frame(
    mz = 100.00005,
    reference = 0,
    sample = 0,
    check.names = FALSE
  )
  composite$`_mz_landmarks_` <- integer()

  CompositeMap_build_composite_tracks(composite)

  expect_length(composite$composite_mass_tracks, 1L)
  expect_equal(
    composite$composite_mass_tracks[[1L]]$intensity,
    c(0, 15, 35, 15, 0)
  )
})

# Verify sample RT reverse mapping, 0-based MassGrid row, and feature table sample columns.
test_that("feature extraction and table generation preserve Python indexing", {
  sample <- make_constructor_sample(
    "sample", 100, list(c(0, 10, 20, 30, 20, 10, 0)), TRUE
  )
  sample$reverse_rt_cal_dict <- c(`1` = 2, `5` = 4)
  experiment <- make_constructor_experiment(list(sample = sample))
  experiment$all_samples <- list(sample)
  composite <- CompositeMap__init__(experiment)
  composite$MassGrid <- data.frame(mz = 100, sample = 0, check.names = FALSE)
  composite$FeatureList <- list(list(
    id_number = "F0",
    parent_masstrack_id = 0,
    mz = 100,
    apex = 3,
    left_base = 1,
    right_base = 5
  ))

  areas <- CompositeMap_extract_features_per_sample(
    composite, sample, CompositeMap_get_peak_area_sum
  )
  CompositeMap_generate_feature_table(composite)

  # The reverse mapping converts [1,5] into a Python closed interval [2,4] with an area of 20+30+20.
  expect_equal(areas, 70)
  expect_equal(composite$FeatureTable$sample, 70)
  expect_equal(composite$FeatureTable$id_number, "F0")
})

# After verifying that the number of landmarks strictly exceeds the threshold, the calibration function is called and the sample status is written back.
test_that("sample RT calibration selects peaks and writes both mappings", {
  scans <- 0:30
  peak_track <- as.integer(
    1000 + 200000 * exp(-((scans - 15)^2) / (2 * 4^2))
  )
  reference <- make_constructor_sample(
    "reference", c(100, 101), list(peak_track, peak_track), TRUE
  )
  sample <- make_constructor_sample(
    "sample", c(100, 101), list(peak_track, peak_track), FALSE
  )
  experiment <- make_constructor_experiment(
    list(reference = reference, sample = sample)
  )
  composite <- CompositeMap(experiment)
  composite$MassGrid <- data.frame(
    mz = c(100, 101),
    reference = c(0, 1),
    sample = c(0, 1),
    check.names = FALSE
  )
  composite$good_reference_landmark_peaks <- list(
    list(apex = 15L, ref_id_num = 0L),
    list(apex = 15L, ref_id_num = 1L)
  )
  fake_calibration <- function(...) {
    list(c(`10` = 11), c(`11` = 10))
  }

  CompositeMap_calibrate_sample_RT(
    composite,
    sample,
    sample$list_mass_tracks,
    calibration_fuction = fake_calibration,
    MIN_PEAK_NUM = 0
  )

  expect_true(sample$is_rt_aligned)
  expect_equal(sample$rt_cal_dict, c(`10` = 11))
  expect_equal(sample$reverse_rt_cal_dict, c(`11` = 10))
  expect_equal(sample$rt_landmarks, c(15L, 15L))
})

# Verification reference landmark selectivity, Python track ID 0 exclusion rules and unique peak detection.
test_that("RT reference selection reproduces Python landmark filtering", {
  scans <- 0:30
  peak_track <- as.integer(
    1000 + 200000 * exp(-((scans - 15)^2) / (2 * 4^2))
  )
  reference <- make_constructor_sample(
    "reference",
    c(100, 120, 140, 160, 180),
    rep(list(peak_track), 5),
    TRUE
  )
  experiment <- make_constructor_experiment(list(reference = reference))
  composite <- CompositeMap(experiment)
  composite$MassGrid <- data.frame(
    mz = c(100, 120, 140, 160, 180),
    reference = 0:4,
    check.names = FALSE
  )
  composite$`_mz_landmarks_` <- 1:5

  selected <- CompositeMap_set_RT_reference(composite, 100000)

  # Python's `if ref_ii` skips 0, so only MassGrid rows 1..4 are returned.
  expect_length(selected, 4L)
  expect_equal(vapply(selected, `[[`, numeric(1), "ref_id_num"), 1:4)
  expect_equal(reference$rt_landmarks, rep(15L, 4))
})

# Verifying composite peak detection adds feature ID, RT and generates sample table immediately.
test_that("global peak detection builds experiment features and table", {
  scans <- 0:100
  intensity <- as.integer(
    1000 + 200000 * exp(-((scans - 50)^2) / (2 * 5^2))
  )
  reference <- make_constructor_sample(
    "reference", 100, list(intensity), TRUE
  )
  reference$rt_numbers <- scans
  reference$list_scan_numbers <- scans
  reference$list_retention_time <- scans * 2
  experiment <- make_constructor_experiment(list(reference = reference))
  experiment$number_scans <- 101L
  experiment$all_samples <- list(reference)
  composite <- CompositeMap(experiment)
  composite$MassGrid <- data.frame(
    mz = 100, reference = 0, check.names = FALSE
  )
  composite$composite_mass_tracks <- list(
    `0` = list(id_number = 0, mz = 100, intensity = intensity)
  )

  CompositeMap_global_peak_detection(composite)

  expect_length(composite$FeatureList, 1L)
  expect_equal(composite$FeatureList[[1L]]$id_number, "F0")
  expect_equal(composite$FeatureList[[1L]]$rtime, 100)
  expect_equal(nrow(composite$FeatureTable), 1L)
  expect_true(composite$FeatureTable$reference[[1L]] > 0)
})

# The DIMS path should take the maximum intensity for each sample and each corresponding mass track.
test_that("DIMS feature table uses maximum track intensity", {
  sample <- make_constructor_sample(
    "sample", c(100, 101), list(c(1, 5, 2), c(0, 3, 8)), TRUE
  )
  experiment <- make_constructor_experiment(list(sample = sample))
  experiment$all_samples <- list(sample)
  composite <- CompositeMap__init__(experiment)
  composite$MassGrid <- data.frame(
    mz = c(100, 101), sample = c(0, 1), check.names = FALSE
  )

  CompositeMap_get_DIMS_feature_table(composite)

  expect_equal(composite$FeatureTable$id_number, c("F0", "F1"))
  expect_equal(composite$FeatureTable$sample, c(5, 8))
})
