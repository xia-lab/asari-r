source(file.path("R", "tools_match_features.R"), local = FALSE)

test_that("get_featureList沿用Python的行列编号", {
  infile <- tempfile(fileext = ".tsv")
  writeLines(c("mz\trt", "100\t2", "200\t3"), infile)
  features <- get_featureList(infile)
  expect_identical(vapply(features, `[[`, "", "id"), c("row2", "row3"))
  expect_equal(features[[2L]]$rtime, 3)
})

test_that("一对多和最佳匹配遵守ppm与严格RT容差", {
  list1 <- list(list(id = "a", mz = 100, rtime = 10))
  list2 <- list(
    list(id = "b", mz = 100.0002, rtime = 11),
    list(id = "c", mz = 100.0001, rtime = 12),
    list(id = "edge", mz = 100, rtime = 15)
  )
  expect_identical(list_match_lcms_features(list1, list2, 5, 5)$a, c("b", "c"))
  expect_identical(best_mz_match_lcms_features(list1, list2, 5, 5)$a, "c")
  expect_identical(best_rt_match_lcms_features(list1, list2, 5, 5)$a, "b")
})

test_that("双向最佳匹配和时间单位转换对应Python", {
  one <- list(list(id = "a", mz = 100, rtime = 10))
  two <- list(list(id = "b", mz = 100.0001, rtime = 11))
  result <- bidirectional_best_match(one, two, 5, 5)
  expect_identical(result[[1L]], list(c("a", "b")))
  expect_equal(convert_min2secs(one)[[1L]]$rtime, 600)
  expect_equal(convert_sec2mins(convert_min2secs(one))[[1L]]$rtime, 10)
})

test_that("八个Python def均有R对应", {
  names <- c(
    "get_featureList", "list_match_lcms_features", "bidirectional_match",
    "best_mz_match_lcms_features", "best_rt_match_lcms_features",
    "bidirectional_best_match", "convert_min2secs", "convert_sec2mins"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
