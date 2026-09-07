source(file.path("R", "tools_plot_scan_seq.R"), local = FALSE)

test_that("plot_msn内部强度函数对应Python", {
  peaks <- list(c(50, 100), c(60, 1000))
  expect_equal(.plot_msn_get_median_height(peaks), log10(550))
  expect_equal(.plot_msn_get_max_height(peaks), 1000)
  expect_equal(.plot_msn_get_max_height(list()), 0)
})

test_that("plot_msn写出双面板PDF", {
  ms1 <- list(list(id = "sp1", rtime = 1, peaks = list(c(50, 100))))
  ms2 <- list(list(id = "sp2", rtime = 2, peaks = list(c(40, 50))))
  outfile <- tempfile(fileext = ".pdf")
  plot_msn(ms1, ms2, nspec_plot = 2L, offset_n = 0L, outfile = outfile)
  expect_true(file.exists(outfile))
})

test_that("plot_scan_seq四个Python def含两个内部函数均有对应", {
  names <- c("summarize_ms_file", "plot_msn", ".plot_msn_get_median_height", ".plot_msn_get_max_height")
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
