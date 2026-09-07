test_that("chromatograms全18个Python def包含内部函数都有明确对应", {
  functions <- c(
    "extract_massTracks_", "extract_single_track_fullrt_length",
    "bin_to_mass_tracks", "build_chromatogram_intensity_aware",
    "build_chromatogram_by_mz_clustering", "merge_two_mass_tracks",
    "get_thousandth_bins", "__rough_check_consecutive_scans__",
    "__check_min_peak_height__", "rt_lowess_calibration",
    "clean_rt_calibration_points", "rt_lowess_calibration_debug",
    "__hacked_lowess__", "savitzky_golay_spline", "dwt_rt_calibrate",
    "remap_intensity_track", "smooth_moving_average", "smooth_lowess"
  )
  expect_length(functions, 18L)
  expect_true(all(vapply(functions, function(name) is.function(get(name)), FALSE)))
})

test_that("强度优先色谱构建与Python分组一致", {
  points <- list(c(100.001, 0, 30), c(100.000, 1, 100), c(101, 2, 50))
  groups <- build_chromatogram_intensity_aware(points, 3L, 0.01)
  expect_length(groups, 2L)
  expect_equal(vapply(groups[[1L]], `[[`, 0, 3L), c(100, 30))
})

test_that("extract_single_track_full_rt_length builds one full-length track", {
  bin <- list(
    c(100.0, 0, 10),
    c(100.2, 1, 20),
    c(99.9, 1, 25)
  )

  track <- extract_single_track_full_rt_length(bin, rt_length = 4)

  expect_equal(track[[1]], 99.95)
  expect_equal(track[[2]], c(10, 25, 0, 0))
})

test_that("bin_to_mass_tracks builds one track for a narrow mz range", {
  bin <- list(
    c(100.0002, 1, 20),
    c(100.0000, 0, 10)
  )

  tracks <- bin_to_mass_tracks(bin, rt_length = 3, mz_tolerance_ppm = 5)

  expect_length(tracks, 1)
  expect_equal(tracks[[1]][[1]], 100.00015)
  expect_equal(tracks[[1]][[2]], c(10, 20, 0))
})

test_that("merge_two_mass_tracks averages mz and sums intensities", {
  track_1 <- list(100.0001, c(10, 20, 0, 5))
  track_2 <- list(100.0003, c(0, 15, 30, 5))

  merged <- merge_two_mass_tracks(track_1, track_2)

  expect_equal(merged[[1]], 100.0002)
  expect_equal(merged[[2]], c(10, 35, 30, 10))
})

test_that("get_thousandth_bins merges adjacent bins and filters noise", {
  make_points <- function(mz, scans, intensities) {
    Map(function(scan, intensity) c(mz, scan, intensity), scans, intensities)
  }

  mz_tree <- list(
    `100000` = make_points(100.0002, 10:14, c(200, 400, 1200, 500, 300)),
    `100001` = make_points(100.0012, 10:14, c(300, 500, 1300, 600, 400)),
    `100010` = make_points(100.0102, c(10, 20, 30, 40, 50), rep(2000, 5)),
    `100020` = make_points(100.0202, 10:14, rep(900, 5))
  )

  bins <- get_thousandth_bins(mz_tree)

  expect_length(bins, 1)
  expect_length(bins[[1]], 10)
})

test_that("clean_rt_calibration_points removes duplicates and outliers", {
  pairs <- c(
    lapply(0:19, function(ii) c(ii, ii)),
    list(c(0, 0), c(30, -70))
  )

  cleaned <- clean_rt_calibration_points(pairs)

  expect_equal(cleaned, lapply(0:19, function(ii) c(ii, ii)))
})

test_that("linear interpolation extrapolates both ends", {
  expect_equal(
    linear_interpolate_with_extrapolation(c(1, 2), c(10, 20), 0:3),
    c(0, 10, 20, 30)
  )
})

test_that("rt_lowess_calibration matches expected scan mappings", {
  sample_apexes <- c(10, 30, 50, 70, 90, 110, 130)
  reference_apexes <- c(11, 29, 52, 69, 93, 108, 134)
  calibration <- rt_lowess_calibration(
    lapply(sample_apexes, function(apex) list(apex = apex)),
    lapply(reference_apexes, function(apex) list(apex = apex)),
    0:150,
    0:160,
    1,
    "sample",
    tempdir()
  )

  expect_equal(unname(calibration[[1]]["29"]), 30L)
  expect_equal(unname(calibration[[1]]["96"]), 98L)
  expect_equal(unname(calibration[[2]]["29"]), 28L)
  expect_equal(unname(calibration[[2]]["97"]), 95L)
})

test_that("remap_intensity_track copies and remaps scan intensities", {
  remapped <- remap_intensity_track(
    c(10, 20, 30, 40),
    numeric(6),
    c(`1` = 3, `3` = 5)
  )

  expect_equal(remapped, c(10, 20, 30, 20, 0, 40))
})

test_that("smooth_moving_average matches nearest-edge smoothing", {
  expect_equal(
    smooth_moving_average(1:5, size = 4),
    c(1.25, 1.75, 2.5, 3.5, 4.25)
  )
})

test_that("smooth_lowess reproduces the LOWESS intensity curve", {
  values <- c(0, 1, 1.8, 3.2, 4.1, 5.5, 5.9, 7.2, 8.4, 9)
  expected <- c(
    -0.038765578226, 0.980142253327, 2.009801765114,
    3.080942042368, 4.172514693893, 5.173918857948,
    6.197753482895, 7.185102515519, 8.158376513741,
    9.137881062834
  )

  expect_equal(smooth_lowess(values, frac = 0.6), expected, tolerance = 1e-12)
})

test_that("placeholder RT calibration functions match Python behavior", {
  expect_null(savitzky_golay_spline(list(), list(), integer(), integer()))
  expect_error(
    dwt_rt_calibrate(list(), list(), integer(), integer()),
    "Not implemented"
  )
})
