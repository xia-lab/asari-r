source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_file_io.R"), local = FALSE)
source(file.path("R", "tools_ms2.R"), local = FALSE)

.test_spectra <- list(
  list(ms_level = 1L, rtime = 10, peaks = matrix(c(100, 2000, 200, 500), ncol = 2, byrow = TRUE)),
  list(ms_level = 2L, precursor_mz = 150, rtime = 11, peaks = matrix(c(50, 2000, 149.5, 3000, 60, 1500), ncol = 2, byrow = TRUE)),
  list(ms_level = 3L, precursor_mz = 150, rtime = 12, peaks = matrix(c(40, 2000), ncol = 2))
)

test_that("谱图提取区分MS级别并过滤前体附近峰", {
  result <- extract_all_spectra_form_file(.test_spectra, min_intensity = 1000, MS2_peak_limit = 1L)
  expect_length(result[[1L]], 1L)
  expect_length(result[[2L]], 1L)
  expect_length(result[[3L]], 1L)
  expect_identical(result[[2L]][[1L]]$id, "sp2")
  expect_identical(result[[2L]][[1L]]$peaks, list(c(50, 2000)))
  expect_length(extract_ms2_form_file(.test_spectra, 1000)[[1L]]$peaks, 2L)
})

test_that("谱图筛选和MS1匹配对应Python", {
  spectra <- list(
    list(precursor_mz = 100.0001, rtime = 10, peaks = list(c(50, 10))),
    list(precursor_mz = 100.0002, rtime = 11, peaks = list(c(50, 20), c(60, 30)))
  )
  expect_identical(get_top_n_ms2_spectra(spectra, 1L)[[1L]], spectra[[2L]])
  expect_identical(get_best_ms2_spectrum(spectra), spectra[[2L]])
  expect_identical(get_top_n_peaks(spectra[[2L]]$peaks, 1L), list(c(60, 30)))
  matched <- get_matched_ms2_ms1(list(list(id = "F1", mz = 100, rtime = 10)), spectra, 30, 5)
  expect_identical(matched[[1L]][[1L]], "F1")
  expect_equal(matched[[1L]][[2L]][[2L]], 1L)
})

test_that("数据库匹配可以重组并寻找MS1峰", {
  results <- list(list(list(INCHIKEY = "K", NAME = "C"), "sp1", 100, 10, 0.9, 2L))
  grouped <- regroup_matches_per_file(results)
  expect_named(grouped[[1L]], "K")
  match <- get_ms1_match_to_ms2cluster(
    grouped[[1L]]$K,
    list(list(rtime = 10, peaks = list(c(100.0001, 99)))),
    mz_tol_ppm = 5, rt_tol = 1
  )
  expect_equal(match[[2L]][[2L]], 99)
})

test_that("十个Python def均有R对应", {
  names <- c(
    "extract_all_spectra_form_file", "regroup_matches_per_file",
    "get_ms1_match_to_ms2cluster", "export_table_ms1match_results",
    "extract_ms2_form_file", "get_top_n_ms2_spectra", "get_top_n_peaks",
    "get_best_ms2_spectrum", "get_matched_ms2_ms1", "match_ms2files_to_features"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
