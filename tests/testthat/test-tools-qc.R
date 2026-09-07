source(file.path("R", "tools_ms2.R"), local = FALSE)
source(file.path("R", "gcms.R"), local = FALSE)
source(file.path("R", "tools_qc.R"), local = FALSE)

.qc_spectra <- list(
  list(ms_level = 1L, polarity = 1L, rtime = 1, peaks = matrix(c(100, 20000, 101, 500), ncol = 2, byrow = TRUE)),
  list(ms_level = 2L, polarity = 1L, rtime = 2, peaks = matrix(c(50, 2000), ncol = 2))
)

test_that("QC内部三个函数计算TIC trio和target", {
  expect_identical(.qc_calcTIC(.qc_spectra), list(1, 20500))
  trios <- .qc_extract_trio(.qc_spectra, min_intensity = 1000)
  expect_identical(trios, list(c(20000, 100, 1)))
  found <- .qc_find_targets(trios, list(list("T", 100, 0)))
  expect_length(found$T, 1L)
})

test_that("generate_qc_report写出HTML摘要", {
  outfile <- tempfile(fileext = ".html")
  result <- generate_qc_report(list(.qc_spectra, outfile, list(list("T", 100, 0))))
  expect_identical(result, outfile)
  text <- paste(readLines(outfile), collapse = "")
  expect_match(text, "Total Ion Chromatogram", fixed = TRUE)
  expect_match(text, "T", fixed = TRUE)
})

test_that("qc六个Python def含三个内部函数均有对应", {
  names <- c("asari_qc_plot", "generate_qc_report", ".qc_calcTIC", ".qc_extract_trio", ".qc_find_targets", "get_dataframe_from_file")
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
