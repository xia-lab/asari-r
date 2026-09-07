source(file.path("R", "tools_entropy_search.R"), local = FALSE)

.test_entropy_backend <- function() {
  backend <- new.env(parent = emptyenv())
  backend$entries <- list(
    list(NAME = "A", peaks = matrix(c(50, 1), ncol = 2)),
    list(NAME = "B", peaks = matrix(c(60, 1), ncol = 2))
  )
  backend$clean_spectrum_for_search <- function(precursor_mz, peaks, precursor_ions_removal_da) {
    peaks[peaks[, 1L] < precursor_mz - precursor_ions_removal_da, , drop = FALSE]
  }
  backend$identity_search <- function(precursor_mz, peaks, ms1_tolerance_in_da,
                                      ms2_tolerance_in_da, output_matched_peak_number) {
    list(c(0.2, 0.9), c(1L, nrow(peaks)))
  }
  backend
}

test_that("普通MS2清理和检索保留元数据", {
  backend <- .test_entropy_backend()
  spectra <- list(list(
    id = "sp1", precursor_mz = 100, rtime = 20,
    peaks = list(c(50, 10), c(99, 20))
  ))
  cleaned <- clean_list_ms2spectra(spectra, backend, param)
  expect_length(cleaned, 1L)
  expect_equal(nrow(cleaned[[1L]][[2L]]), 1L)
  matched <- search_ms2_spectra(cleaned, backend, param)
  expect_identical(matched[[1L]][[1L]]$NAME, "B")
  expect_null(matched[[1L]][[1L]]$peaks)
  expect_equal(matched[[1L]][[5L]], 0.9)
})

test_that("MS1分组MS2清理和搜索以feature id返回", {
  backend <- .test_entropy_backend()
  grouped <- list(F1 = list(list(
    precursor_mz = 100,
    peaks = list(c(50, 10), c(99, 20))
  ), 2L))
  cleaned <- get_cleaned_experimental_ms2(grouped, backend, param)
  result <- search_ms1matchedms2_spectra(cleaned, backend, param)
  expect_named(result, "F1")
  expect_identical(result$F1[[1L]]$NAME, "B")
  expect_equal(result$F1[[3L]], 1L)
})

test_that("四个Python def均有R对应", {
  names <- c(
    "clean_list_ms2spectra", "search_ms2_spectra",
    "get_cleaned_experimental_ms2", "search_ms1matchedms2_spectra"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
