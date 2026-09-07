source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_file_io.R"), local = FALSE)
source(file.path("R", "tools_match_features.R"), local = FALSE)
source(file.path("R", "tools_cosine.R"), local = FALSE)
source(file.path("R", "tools_msp_parser.R"), local = FALSE)
source(file.path("R", "gcms.R"), local = FALSE)
source(file.path("R", "annotate.R"), local = FALSE)

test_that("LCMS_Annotation质量校准和孤儿封装对应Python", {
  params <- default_parameters()
  params$outdir <- tempdir()
  annotation <- LCMS_Annotation(params)
  annotation$list_features <- list(
    list(id = "F1", id_number = "F1", mz = 100, peak_area = 10),
    list(id = "F2", id_number = "F2", mz = 200, peak_area = 5)
  )
  database <- new.env(parent = emptyenv())
  database$evaluate_mass_accuracy_ratio <- function(...) 10e-6
  annotation$KCD <- database
  annotation$db_mass_calibrate()
  expect_equal(annotation$list_features[[1L]]$mz, 100 / 1.00001)
  expect_equal(annotation$list_features[[1L]]$mz_corrected_by_division, 1.00001)

  compounds <- annotation$append_orphans_to_epmCpds(list())
  expect_length(compounds, 2L)
  expect_equal(compounds[[1L]]$interim_id, 100000L)
})

test_that("select_unique_compound_features按面积和峰形选择", {
  annotation <- LCMS_Annotation(default_parameters())
  compounds <- list(A = list(
    neutral_formula = "C",
    MS1_pseudo_Spectra = list(
      list(id_number = "F1", peak_area = 10, goodness_fitting = 0.8, ion_relation = "a"),
      list(id_number = "F2", peak_area = 20, goodness_fitting = 0.7, ion_relation = "b")
    )
  ))
  annotation$select_unique_compound_features(compounds)
  expect_named(annotation$selected_unique_features, "F2")
  expect_identical(annotation$selected_unique_features$F2[[3L]], "b")
})

test_that("标准库匹配可压缩并导出表和JSON", {
  features <- list(F1 = list(id = "F1", mz = 100, rtime = 10))
  libraries <- list(
    v2r2024_A = list(name = "Auth", mz = 100, rtime = 10, identifier = "I", ion = "M+H", isotope = "M0"),
    r1_B = list(ion_csm = "M+H", top_recommendation_name = "Pred", top_recommendation_score = 0.8, HMDB = "H")
  )
  matched <- list(F1 = c("v2r2024_A", "r1_B"))
  concise <- get_concise_annotation("F1", matched$F1, features, libraries)
  expect_identical(concise$lib_name, "Auth")
  expect_identical(concise$CSM_top_recommendation_name, "Pred")
  outfile <- tempfile(fileext = ".tsv")
  export_combined_anno_table(matched, features, libraries, outfile)
  expect_equal(length(readLines(outfile)), 2L)
})

test_that("annotate全部15个Python def含重复定义均有R对应", {
  names <- c(
    "annotate_project", "annotate_gcms_full", "LCMS_Annotation__init__",
    "LCMS_Annotation_annotate_user_featuretable", "LCMS_Annotation_load_annotation_db",
    "LCMS_Annotation_db_mass_calibrate_documented", "LCMS_Annotation_db_mass_calibrate",
    "LCMS_Annotation_append_orphans_to_epmCpds", "LCMS_Annotation_export_peak_annotation",
    "LCMS_Annotation_select_unique_compound_features", "annotate_user_featuretable",
    "annoate_by_standards", "get_concise_annotation",
    "export_combined_anno_table", "export_combined_anno_json"
  )
  expect_length(names, 15L)
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
