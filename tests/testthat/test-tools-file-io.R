source(file.path("R", "tools_file_io.R"), local = FALSE)

test_that("read_features_from_asari_table保留字段和值", {
  header <- paste(c(paste0("c", 0:10), "sampleA", "sampleB"), collapse = "\t")
  row <- paste(c("F1", "100.25", "30", "28", "33", "M1", "5000", "0.9", "0.8", "12", "2", "10", "20"), collapse = "\t")
  result <- suppressWarnings(read_features_from_asari_table(paste(header, row, sep = "\n")))
  expect_equal(result[[1L]], 2L)
  expect_length(result[[2L]], 1L)
  expect_identical(result[[2L]][[1L]]$id, "F1")
  expect_equal(result[[2L]][[1L]]$mz, 100.25)
  expect_equal(result[[2L]][[1L]]$apex, 30)
  expect_equal(result[[2L]][[1L]]$detection_counts, 2L)
})

test_that("JSON风格记录可以导出再按字符串读回", {
  records <- list(
    list(name = "alpha", id = "x1", score = 2),
    list(name = "beta", id = "x2", score = 3)
  )
  outfile <- tempfile(fileext = ".tsv")
  export_json_to_table(records, outfile)
  expect_identical(readLines(outfile)[[1L]], "id\tname\tscore")
  restored <- read_table_to_json(outfile)
  expect_identical(restored[[1L]], list(id = "x1", name = "alpha", score = "2"))
})

test_that("三个Python def均有R对应", {
  expect_true(is.function(read_features_from_asari_table))
  expect_true(is.function(export_json_to_table))
  expect_true(is.function(read_table_to_json))
})
