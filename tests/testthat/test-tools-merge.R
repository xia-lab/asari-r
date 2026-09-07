source(file.path("R", "tools_merge.R"), local = FALSE)

test_that("bin_by_median按动态中位数分箱", {
  tuples <- list(list(1, "a"), list(1.1, "b"), list(3, "c"))
  # Python returns list-of-lists. R also retains each bin as a list and does not automatically simplify it to a character vector.
  expect_identical(
    bin_by_median(tuples, function(value) 0.5),
    list(list("a", "b"), list("c"))
  )
  expect_identical(`_formula_selectivity_`(c("A", "A", "B")), list(A = 2L, B = 1L))
})

test_that("读取表格并在同一质量轨迹内对齐RT", {
  infile <- tempfile(fileext = ".tsv")
  writeLines(c(
    "feature_id\tformula_mass\tmz\trtime\tarea",
    "100@10\tM100\t100\t10\t1",
    "100@30\tM100\t100\t30\t2"
  ), infile)
  parsed <- `_read_asari_ftables_`(infile)
  expect_equal(parsed[[1L]][["100@10"]][[2L]], 100)
  grouped <- `_masstrace2features_`(list(c("100@10", "100@30"), c("100@11")), 15)
  expect_length(grouped, 2L)
  expect_identical(grouped[[1L]], list("100@10", "100@11"))
})

test_that("merge_feature_tables生成包含来源信息的结果", {
  first <- tempfile(fileext = ".tsv")
  second <- tempfile(fileext = ".tsv")
  outfile <- tempfile(fileext = ".tsv")
  writeLines(c("feature_id\tformula_mass\tmz\trtime", "100@10\tM100\t100\t10"), first)
  writeLines(c("feature_id\tformula_mass\tmz\trtime", "100@11\tM100\t100\t11"), second)
  merge_feature_tables(c(first, second), outfile = outfile)
  lines <- readLines(outfile)
  expect_match(lines[[1L]], "new_feature_id", fixed = TRUE)
  expect_match(lines[[2L]], "100@10,100@11", fixed = TRUE)
})

test_that("源码中的旧接口和实际def均有对应", {
  names <- c(
    "match2", "__bin_by_median_rt__", "_formula_selectivity_",
    "tools_merge_bin_by_median", "_read_asari_ftables_", "_masstrace2features_",
    ".merge_find_min_delta", "merge_feature_tables"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
