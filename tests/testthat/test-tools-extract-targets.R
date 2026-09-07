source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_extract_targets.R"), local = FALSE)

test_that("export_targets按ppm匹配并在首列写入查询质量", {
  outfile <- tempfile(fileext = ".tsv")
  parameters <- list(
    target = c(100.0002, 300), mz_tolerance_ppm = 5,
    outdir = tempdir(), output_feature_table = "Feature_table.tsv"
  )
  table <- data.frame(id = c("F1", "F2"), mz = c(100, 200), check.names = FALSE)
  export_targets(parameters, table, 1L, outfile)
  result <- utils::read.delim(outfile, check.names = FALSE)
  expect_equal(nrow(result), 1L)
  expect_equal(result$query_target, 100.0002)
  expect_identical(result$id, "F1")
})

test_that("空target与Python一样不导出", {
  outfile <- tempfile()
  export_targets(list(target = NULL), data.frame(mz = 1), 1L, outfile)
  expect_false(file.exists(outfile))
})

test_that("唯一Python def有R对应", {
  expect_true(is.function(export_targets))
})
