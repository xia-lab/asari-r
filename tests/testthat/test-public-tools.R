# The public tool test only loads the modules required by this file to ensure that each entry can be independently audited.
source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "asari_process.R"), local = FALSE)
source(file.path("R", "analyze.R"), local = FALSE)
source(file.path("R", "tools_match_features.R"), local = FALSE)
source(file.path("R", "tools_ms2.R"), local = FALSE)
source(file.path("R", "tools_qc.R"), local = FALSE)
source(file.path("R", "tools_plot_scan_seq.R"), local = FALSE)
source(file.path("R", "public_tools.R"), local = FALSE)

test_that("DIMS公开入口固定工作流并关闭RT对齐", {
  original <- asari_process
  on.exit(assign("asari_process", original, envir = .GlobalEnv), add = TRUE)
  assign("asari_process", function(...) list(...), envir = .GlobalEnv)

  value <- asari_process_dims("sample.mzML", parameters = list(min_timepoints = 3L))
  expect_identical(value$rt_align, FALSE)
  expect_identical(value$parameters$workflow, "DIMS")
  expect_identical(value$parameters$rt_align_on, FALSE)
  expect_identical(value$parameters$min_timepoints, 3L)
})

test_that("单文件分析返回可继续使用的质量轨迹表", {
  skip_if_not_installed("mzR")
  infile <- tempfile(fileext = ".mzML")
  file.create(infile)
  header <- data.frame(
    msLevel = c(1L, 1L, 2L), polarity = c(1L, 1L, 1L),
    retentionTime = c(1, 2, 2)
  )
  extractor <- function(...) list(
    tracks = list(list(100, c(0, 10)), list(101.003355, c(0, 20))),
    rt_numbers = c(0L, 1L), rt_times = c(1, 2)
  )
  pair_finder <- function(...) list(c(0L, 1L))
  value <- asari_analyze(
    infile, min_intensity = 1, min_timepoints = 1, min_peak_height = 1,
    parameters = list(
      scan_header = header,
      extract_mass_tracks = extractor,
      find_mzdiff_pairs = pair_finder
    )
  )

  expect_s3_class(value, "asari_analysis")
  expect_identical(value$ms1_scans, 2L)
  expect_identical(value$ms2_scans, 1L)
  expect_identical(value$isotope_pair_count, 1L)
  expect_equal(value$mass_track_table$max_intensity, c(10, 20))
})

test_that("项目结果读取不要求dashboard pickle", {
  project <- tempfile("asari-project-")
  dir.create(file.path(project, "export"), recursive = TRUE)
  jsonlite::write_json(
    list(number_of_samples = 1L, output_feature_table = "Custom.tsv"),
    file.path(project, "project.json")
  )
  table <- data.frame(
    id_number = "F0", mz = 100, rtime = 10,
    rtime_left_base = 9, rtime_right_base = 11,
    parent_masstrack_id = 0, peak_area = 1000, cSelectivity = 1,
    goodness_fitting = 1, snr = 10, detection_counts = 1,
    sample = 1000, check.names = FALSE
  )
  utils::write.table(
    table, file.path(project, "preferred_Custom.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    table, file.path(project, "export", "full_Custom.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  value <- asari_read_results(project)
  expect_s3_class(value, "asari_project_results")
  expect_equal(nrow(value$preferred), 1L)
  expect_equal(nrow(value$full), 1L)
})

test_that("公开表格工具完成双向匹配和目标提取", {
  left <- data.frame(
    id_number = c("L1", "L2"), mz = c(100, 200), rtime = c(10, 20)
  )
  right <- data.frame(
    id_number = c("R1", "R2"), mz = c(100.0002, 300), rtime = c(11, 30)
  )
  matched <- asari_compare_features(left, right, ppm = 5, rt_tolerance = 5)
  expect_equal(nrow(matched), 1L)
  expect_identical(matched$left_id, "L1")
  expect_identical(matched$right_id, "R1")

  targeted <- asari_extract_targets(left, c(100.0002, 400), ppm = 5)
  expect_equal(nrow(targeted), 1L)
  expect_identical(targeted$id_number, "L1")
  expect_equal(attr(targeted, "unmatched_targets"), 400)
})

test_that("QC和扫描摘要公开入口写出明确目标文件", {
  input <- tempfile(fileext = ".mzML")
  file.create(input)
  feature_pdf <- tempfile(fileext = ".pdf")
  feature_table <- data.frame(
    snr = c(4, 8), peak_area = c(1000, 2000),
    cSelectivity = c(0.8, 1), goodness_fitting = c(0.9, 1)
  )
  expect_true(file.exists(asari_feature_qc(feature_table, feature_pdf)))

  original <- .ms2_read_spectra
  on.exit(assign(".ms2_read_spectra", original, envir = .GlobalEnv), add = TRUE)
  spectra <- list(
    list(
      id = "sp1", ms_level = 1L, polarity = 1L, rtime = 1,
      peaks = matrix(c(100, 20000, 101, 1000), ncol = 2, byrow = TRUE)
    ),
    list(
      id = "sp2", ms_level = 2L, polarity = 1L, precursor_mz = 100,
      rtime = 2, peaks = matrix(c(50, 2000), ncol = 2)
    )
  )
  assign(".ms2_read_spectra", function(infile) spectra, envir = .GlobalEnv)
  html <- tempfile(fileext = ".html")
  pdf <- tempfile(fileext = ".pdf")
  expect_true(file.exists(asari_qc_report(input, html, spikeins = list(list("T", 100, 0)))))
  expect_true(file.exists(asari_scan_summary(input, pdf, nspec_plot = 2L, offset_n = 0L)))
})

test_that("工作流状态列出新增公开入口和必要输入", {
  status <- asari_available_workflows()
  expect_true(status$direct_use[status$workflow == "LC"])
  expect_true(status$direct_use[status$workflow == "DIMS"])
  expect_true(status$direct_use[status$workflow == "GC"])
  expect_match(status$note[status$workflow == "GC"], "Kovats", fixed = TRUE)
  expect_true(status$direct_use[status$workflow == "LCMSMS"])
})
