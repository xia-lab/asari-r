source(file.path("R", "tools_msp_parser.R"), local = FALSE)

test_that("简单MSP解析和字段标准化对应Python", {
  file <- tempfile(fileext = ".msp")
  writeLines(c("Name: Test", "MW: 100", "Num Peaks: 2", "50 10; 60 20"), file)
  parsed <- parse_msp_to_listdict(file)
  expect_identical(parsed[[1L]]$Name, "Test")
  expect_length(parsed[[1L]]$peaks, 2L)
  standardized <- msp_standarize(parsed, MSP_dict)
  expect_identical(standardized[[1L]]$ExactMass, "100")
})

test_that("MGF解析保留PEPMASS CHARGE和峰数组", {
  file <- tempfile(fileext = ".mgf")
  writeLines(c(
    "BEGIN IONS", "TITLE=test", "PEPMASS=100.5 200", "CHARGE=2+",
    "50 10", "bad peak", "END IONS"
  ), file)
  spectra <- parse_mgf(file)
  expect_length(spectra, 1L)
  expect_identical(spectra[[1L]]$params$PEPMASS, list(100.5, 200))
  expect_equal(spectra[[1L]]$params$CHARGE, 2L)
  expect_equal(spectra[[1L]]$mzs, 50)
})

test_that("宽松MSP峰解析兼容空格逗号和冒号", {
  expect_identical(
    parse_peak_token("100:45:loss:H2O"),
    list(mz = 100, intensity = 45, annotation = "loss:H2O")
  )
  peaks <- parse_peak_line("100 45, 120:30 140 10 note")
  # The original Python parser will regard 120:30 immediately following the space pair as an annotation of the previous peak.
  expect_equal(vapply(peaks, `[[`, 0, "mz"), c(100, 140))
  expect_identical(peaks[[1L]]$annotation, "120:30")
  expect_identical(peaks[[2L]]$annotation, "note")

  file <- tempfile(fileext = ".msp")
  writeLines(c("Name: A", "Num Peaks: 2", "100 1", "120:2"), file)
  parsed <- parse_msp(file)
  expect_length(parsed, 1L)
  expect_length(parsed[[1L]]$peaks, 2L)
})

test_that("八个Python def含两个内部函数均有R对应", {
  names <- c(
    "msp_standarize", "parse_msp_to_listdict", "parse_mgf",
    ".parse_mgf_pepmass", ".parse_mgf_charge", "parse_peak_token",
    "parse_peak_line", "parse_msp"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
