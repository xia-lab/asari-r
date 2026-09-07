source(file.path("R", "package_init.R"), local = FALSE)

test_that("包初始化版本和Python一致", {
  expect_identical(`__version__`, "1.17.0")
  expect_identical(asari_version, `__version__`)
  expect_length(ASARI_CORE_MODULES, 12L)
  expect_true(all(c("main", "workflow", "peaks") %in% ASARI_CORE_MODULES))
})
