source(file.path("R", "db_init.R"), local = FALSE)

test_that("数据库兼容层完整覆盖三个Python资源", {
  expect_setequal(
    ASARI_DB_RESOURCES,
    c("gcms_libraries.json", "mass_indexed_compounds.pickle", "emp_cpds_trees.pickle")
  )
  path <- asari_db_path("gcms_libraries.json")
  expect_true(file.exists(path))
  expect_match(path, "gcms_libraries.json", fixed = TRUE)
})

test_that("JSON数据库资源可以读取", {
  skip_if_not_installed("jsonlite")
  library <- load_asari_db_resource("gcms_libraries.json")
  expect_true(is.list(library))
  expect_gt(length(library), 0L)
})
