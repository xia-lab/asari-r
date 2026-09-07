source(file.path("R", "default_parameters.R"), local = FALSE)

test_that("default_parameters完整对应Python默认配置", {
  params <- default_parameters()

  expect_equal(length(params) - 1L, 55L)
  expect_identical(params$project_name, "asari_project")
  expect_identical(params$database_mode, "ondisk")
  expect_equal(params$mass_range, c(50, 2000))
  expect_equal(params$mz_tolerance_ppm, 5)
  expect_equal(params$min_timepoints, 6L)
  expect_equal(params$min_intensity_threshold, 1000)
  expect_equal(params$min_peak_height, 100000)
  expect_identical(params$workflow, "LC")
  expect_true(params$rt_align_on)
  expect_false(params$anno)
  expect_identical(params$output_feature_table, "Feature_table.tsv")

  # Compatible fields do not change the value of the Python parameter itself.
  expect_equal(params$min_intensity, params$min_intensity_threshold)
})

test_that("每次调用返回独立参数列表", {
  first <- default_parameters()
  second <- default_parameters()
  first$mode <- "neg"
  expect_identical(second$mode, "pos")
  expect_identical(PARAMETERS$mode, "pos")
})

test_that("加合物和同位素规则完整保留", {
  expect_equal(PROTON, 1.00727646677)
  expect_equal(electron, 0.000549)
  expect_length(adduct_search_patterns_pos, 4L)
  expect_identical(adduct_search_patterns, adduct_search_patterns_pos)
  expect_length(adduct_search_patterns_neg, 6L)
  expect_length(isotope_search_patterns, 4L)
  expect_length(extended_adducts, 17L)
  expect_equal(isotope_search_patterns[[1L]][[1L]], 1.003355)
  expect_match(readme_doc_str, "preferred_Feature_table.tsv", fixed = TRUE)
})
