source(file.path("R", "tools_cosine.R"), local = FALSE)

test_that("find_matches返回Python风格0基索引", {
  matches <- find_matches(c(50, 80, 100), c(49.99, 80.01, 120), 0.05)
  expect_identical(matches[[1L]], c(0L, 1L))
  expect_identical(matches[[2L]], c(0L, 1L))
})

test_that("collect_peak_pairs支持偏移搜索和最小匹配数", {
  ref <- matrix(c(50, 1, 80, 0.5), ncol = 2, byrow = TRUE)
  qry <- matrix(c(55, 1, 85, 0.5), ncol = 2, byrow = TRUE)
  none <- collect_peak_pairs(ref, qry, 1L, TRUE, 0.05, 0)
  shifted <- collect_peak_pairs(ref, qry, 2L, TRUE, 0.05, 5)
  expect_length(none[[1L]], 0L)
  expect_identical(shifted[[1L]], c(0L, 1L))
  expect_identical(shifted[[2L]], c(0L, 1L))
})

test_that("cosine_similarity复现标准和反向余弦行为", {
  ref <- matrix(c(50, 8, 70, 100, 80, 50, 100, 50), ncol = 2, byrow = TRUE)
  qry <- matrix(c(55, 38, 80, 66, 90, 999), ncol = 2, byrow = TRUE)
  standard <- cosine_similarity(qry, ref, tolerance = 0.05, penalty = 0)
  reverse <- cosine_similarity(qry, ref, tolerance = 0.05, penalty = 1)
  expect_equal(standard[[2L]], 1L)
  expect_equal(reverse[[2L]], 1L)
  expect_gt(reverse[[1L]], standard[[1L]])
  expect_lte(reverse[[1L]], 1)
})

test_that("四个Python def均有R对应", {
  expect_true(all(vapply(
    c("find_matches", "collect_peak_pairs", "score_matches", "cosine_similarity"),
    function(name) is.function(get(name)), FALSE
  )))
})
