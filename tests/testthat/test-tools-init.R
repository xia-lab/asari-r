source(file.path("R", "tools_init.R"), local = FALSE)

test_that("空tools初始化模块不额外定义Python def", {
  expect_length(ASARI_TOOL_MODULES, 16L)
  expect_true(all(c("cosine", "gc_annotation", "qc") %in% ASARI_TOOL_MODULES))
})
