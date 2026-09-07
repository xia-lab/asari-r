source(file.path("R", "tools_ms2.R"), local = FALSE)
source(file.path("R", "tools_plot.R"), local = FALSE)

test_that("get_plot_region_from_file使用0基扫描和严格范围", {
  spectra <- list(
    list(ms_level = 1L, peaks = matrix(c(100, 10), ncol = 2)),
    list(ms_level = 1L, peaks = matrix(c(100, 10, 110, 20), ncol = 2, byrow = TRUE)),
    list(ms_level = 2L, peaks = matrix(c(100, 30), ncol = 2))
  )
  result <- get_plot_region_from_file(spectra, 0, 2, 99, 105, 1L)
  expect_identical(result, list(c(1, 100, 10)))
})

test_that("mirror_plot和选择性图写出PDF", {
  mirror <- tempfile(fileext = ".pdf")
  mirror_plot(matrix(c(50, 10), ncol = 2), matrix(c(50.001, 20), ncol = 2), match_tol = 0.01, outfile = mirror)
  expect_true(file.exists(mirror))
  prefix <- tempfile()
  plot_mSelectivity(c(50, 60), c(1, 0.5), outfile = prefix)
  expect_true(file.exists(paste0(prefix, ".pdf")))
})

test_that("plot全部13个Python def均有R对应", {
  names <- c(
    "get_plot_region_from_file", "plot_scatter_map_region",
    "double_scatter_map_region", "with_line_scatter_map_region",
    "plot_masstrack", "mirror_plot", "plot_mSelectivity", "plot_cSelectivity",
    "pca_ftable", "plot_correlations", "plot_peaks_masstrace",
    "plot_peaks", "plot_sample_rt_calibration"
  )
  expect_length(names, 13L)
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
