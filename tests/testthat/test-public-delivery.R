# Public delivery tests use small explicit inputs and do not connect to the network or hide data.
source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_file_io.R"), local = FALSE)
source(file.path("R", "tools_match_features.R"), local = FALSE)
source(file.path("R", "tools_cosine.R"), local = FALSE)
source(file.path("R", "tools_msp_parser.R"), local = FALSE)
source(file.path("R", "tools_ms2.R"), local = FALSE)
source(file.path("R", "gcms.R"), local = FALSE)
source(file.path("R", "public_tools.R"), local = FALSE)
source(file.path("R", "public_annotation_ms2.R"), local = FALSE)
source(file.path("R", "public_interfaces.R"), local = FALSE)

test_that("LC公开注释同时支持离子m/z和中性质量", {
  features <- data.frame(
    id_number = c("F1", "F2"), mz = c(100.0002, 101.00727646677),
    rtime = c(10, 20)
  )
  database <- data.frame(
    id = c("STD", "NEUTRAL"), name = c("Standard", "Neutral"),
    mz = c(100, NA), neutral_mass = c(NA, 100), rtime = c(10, NA)
  )
  outfile <- tempfile(fileext = ".tsv")
  result <- asari_annotate_lc(
    features, database, outfile, mode = "pos", ppm = 5,
    rt_tolerance = 2, include_adducts = FALSE
  )
  expect_s3_class(result, "asari_annotation")
  expect_setequal(result$compound_id, c("STD", "NEUTRAL"))
  expect_true(file.exists(outfile))
})

test_that("GC公开入口校验路径并返回固定输出名", {
  files <- replicate(3L, tempfile(), simplify = TRUE)
  invisible(vapply(files, file.create, FALSE))
  outdir <- tempfile("gc-annotation-")
  original <- annotate_gcms_full
  on.exit(assign("annotate_gcms_full", original, envir = .GlobalEnv), add = TRUE)
  assign("annotate_gcms_full", function(...) invisible(NULL), envir = .GlobalEnv)
  result <- asari_annotate_gc(files[[1L]], files[[2L]], files[[3L]], outdir, "demo")
  expect_s3_class(result, "asari_gc_annotation")
  expect_match(result$feature_annotations, "Features_demo.tsv", fixed = TRUE)
})

test_that("MS2提取、MS1匹配和本地谱库检索有公开结果", {
  infile <- tempfile(fileext = ".mzML")
  file.create(infile)
  spectra <- list(
    list(ms_level = 1L, rtime = 10, peaks = matrix(c(100, 5000), ncol = 2)),
    list(
      id = "sp2", ms_level = 2L, precursor_mz = 100, rtime = 10,
      peaks = matrix(c(40, 4000, 60, 3000, 99.5, 5000), ncol = 2, byrow = TRUE)
    )
  )
  original <- .ms2_read_spectra
  on.exit(assign(".ms2_read_spectra", original, envir = .GlobalEnv), add = TRUE)
  assign(".ms2_read_spectra", function(path) spectra, envir = .GlobalEnv)
  extracted <- asari_extract_ms2(infile, min_intensity = 1000)
  expect_length(extracted$ms2, 1L)
  library <- list(list(
    id = "L1", name = "Library one", precursor_mz = 100,
    peaks = list(c(40, 100), c(60, 80))
  ))
  searched <- asari_search_ms2(extracted, library, min_similarity = 0, method = "cosine")
  expect_equal(nrow(searched), 1L)
  expect_identical(searched$library_id, "L1")
  expect_true(searched$matched_peaks >= 2L)
})

test_that("离线Dashboard只需项目表格并生成可搜索HTML", {
  project <- tempfile("asari-dashboard-")
  dir.create(file.path(project, "export"), recursive = TRUE)
  jsonlite::write_json(
    list(project_name = "demo", number_of_samples = 1L, output_feature_table = "Feature_table.tsv"),
    file.path(project, "project.json"), auto_unbox = TRUE
  )
  table <- data.frame(
    id_number = "F0", mz = 100, rtime = 10, rtime_left_base = 9,
    rtime_right_base = 11, parent_masstrack_id = 1, peak_area = 1000,
    cSelectivity = 1, goodness_fitting = 1, snr = 10,
    detection_counts = 1, sample = 1000, check.names = FALSE
  )
  utils::write.table(table, file.path(project, "preferred_Feature_table.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(table, file.path(project, "export", "full_Feature_table.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  outfile <- asari_dashboard(project, output = tempfile(fileext = ".html"))
  html <- paste(readLines(outfile), collapse = "\n")
  expect_match(html, "Feature map", fixed = TRUE)
  expect_match(html, "Search any visible value", fixed = TRUE)
})

test_that("GUI向导可先返回配置而不启动计算", {
  configuration <- asari_gui(
    input = "sample.mzML", output = "results", workflow = "DIMS", run = FALSE
  )
  expect_identical(configuration$workflow, "DIMS")
  expect_match(configuration$input, "sample.mzML", fixed = TRUE)
})

test_that("CLI列出工作流时不要求输入文件", {
  output <- capture.output(result <- asari_cli("list_workflows"))
  expect_true(any(grepl("asari_process", output, fixed = TRUE)))
  expect_true(any(result$workflow == "dashboard"))
})

test_that("RAW公开入口可使用注入运行器而不连网测试", {
  raw <- tempfile(fileext = ".raw")
  file.create(raw)
  original <- mzMLconverter
  on.exit(assign("mzMLconverter", original, envir = .GlobalEnv), add = TRUE)
  assign("mzMLconverter", function(...) {
    list(command_template = c("fake-converter", "-i", "INPUT", "-o", "OUTPUT"))
  }, envir = .GlobalEnv)
  old_runner <- getOption("asariR.converter_run")
  on.exit(options(asariR.converter_run = old_runner), add = TRUE)
  options(asariR.converter_run = function(command) 0L)
  result <- asari_convert_raw(raw, output_dir = tempfile("raw-output-"), multicores = 1L)
  expect_false(result$dry_run)
  expect_match(result$output_files, ".mzML", fixed = TRUE)
})
