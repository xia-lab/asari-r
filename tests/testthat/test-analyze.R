source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "analyze.R"), local = FALSE)

test_that("match_mzdiff_pairs_by_rt按0基apex和严格窗口筛选", {
  tracks <- list(
    list(id_number = 0L, intensity = c(0, 10, 0)),
    list(id_number = 1L, intensity = c(0, 0, 10)),
    list(id_number = 2L, intensity = c(rep(0, 20), 10))
  )
  pairs <- list(c(0L, 1L), c(0L, 2L))
  expect_identical(match_mzdiff_pairs_by_rt(pairs, tracks, 100, min_scans_window = 10), list(c(0L, 1L)))
})

test_that("get_file_masstrack_stats使用注入提取器并返回landmark", {
  parameters <- list(
    mz_tolerance_ppm = 5, min_intensity_threshold = 100,
    min_timepoints = 5L, min_peak_height = 1000,
    scan_header = data.frame(msLevel = c(1L, 2L), polarity = c(1L, NA)),
    extract_mass_tracks = function(...) list(
      rt_numbers = 0:2, rt_times = c(0, 1, 2),
      tracks = list(list(100, c(0, 10, 0)), list(101.003355, c(0, 5, 0)))
    ),
    find_mzdiff_pairs = function(...) list(c(0L, 1L))
  )
  result <- get_file_masstrack_stats("x.mzML", parameters)
  expect_identical(result[[1L]], c(0L, 1L))
  expect_identical(result[[2L]], "pos")
  expect_equal(result[[3L]], 5L)
  sample <- get_file_masstrack_stats("x.mzML", parameters, TRUE)
  expect_length(sample$list_mass_tracks, 2L)
})

test_that("estimate_min_peak_height聚合成功结果并实施下限", {
  old <- get("get_file_masstrack_stats", envir = .GlobalEnv)
  assign("get_file_masstrack_stats", function(infile, parameters, return_sample = FALSE) {
    list(integer(), "pos", if (infile == "a") 100L else 200L)
  }, envir = .GlobalEnv)
  on.exit(assign("get_file_masstrack_stats", old, envir = .GlobalEnv), add = TRUE)
  parameters <- list(
    min_peak_height = 1000, min_min_peak_height = 10,
    dynamic_range = 1000, num_files_to_check = NULL, multicores = 1L
  )
  expect_equal(estimate_min_peak_height(c("a", "b"), parameters), 75L)
})

test_that("analyze五个Python def均有R对应", {
  names <- c(
    "analyze_single_sample", "get_file_masstrack_stats",
    "match_mzdiff_pairs_by_rt", "__wrapped_get_file_masstrack_stats",
    "estimate_min_peak_height"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
