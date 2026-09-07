# Merge, feature map, and plot public entries only use explicit temporary files.
source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_ms2.R"), local = FALSE)
source(file.path("R", "tools_plot.R"), local = FALSE)
source(file.path("R", "tools_merge.R"), local = FALSE)
source(file.path("R", "tools_feature_graph.R"), local = FALSE)
source(file.path("R", "public_tools.R"), local = FALSE)
source(file.path("R", "public_analysis.R"), local = FALSE)

.public_analysis_table <- function(ids, mz, rt, sample_name, values) {
  result <- data.frame(
    id_number = ids, mz = mz, rtime = rt,
    rtime_left_base = rt - 1, rtime_right_base = rt + 1,
    parent_masstrack_id = seq_along(ids), peak_area = values,
    cSelectivity = 1, goodness_fitting = 1, snr = 10,
    detection_counts = 1, check.names = FALSE
  )
  result[[sample_name]] <- values
  result
}

test_that("标准特征表合并保留来源ID和样本强度", {
  left <- .public_analysis_table(c("F1", "F2"), c(100, 200), c(10, 20), "S1", c(10, 20))
  right <- .public_analysis_table(c("G1", "G2"), c(100.0002, 300), c(11, 30), "S2", c(11, 30))
  paths <- c(tempfile(fileext = ".tsv"), tempfile(fileext = ".tsv"))
  utils::write.table(left, paths[[1L]], sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(right, paths[[2L]], sep = "\t", row.names = FALSE, quote = FALSE)
  result <- asari_merge_feature_tables(
    paths, tempfile(fileext = ".tsv"), ppm = 5, rt_tolerance = 5
  )
  expect_equal(nrow(result), 3L)
  expect_true(any(result$detection_tables == 2L))
  expect_true(any(grepl("F1", result$source_ids, fixed = TRUE)))
  expect_true(any(grepl("::S1", names(result), fixed = TRUE)))
})

test_that("GC特征图返回共洗脱谱簇", {
  table <- .public_analysis_table(
    c("F1", "F2", "F3"), c(50, 60, 70), c(10, 10.2, 30),
    "S1", c(100, 80, 40)
  )
  path <- tempfile(fileext = ".tsv")
  utils::write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE)
  result <- asari_feature_graph(path, rt_tolerance = 0.5)
  expect_setequal(result$feature_id, c("F1", "F2"))
  expect_equal(length(unique(result$cluster_id)), 1L)
})

test_that("镜像图和mzML区域图写出PDF", {
  mirror <- asari_mirror_plot(
    matrix(c(50, 100, 60, 50), ncol = 2, byrow = TRUE),
    matrix(c(50, 80, 70, 20), ncol = 2, byrow = TRUE),
    tempfile(fileext = ".pdf"), tolerance = 0.01
  )
  expect_true(file.exists(mirror))

  infile <- tempfile(fileext = ".mzML")
  file.create(infile)
  original <- .ms2_read_spectra
  on.exit(assign(".ms2_read_spectra", original, envir = .GlobalEnv), add = TRUE)
  assign(".ms2_read_spectra", function(path) list(
    list(ms_level = 1L, peaks = matrix(c(100, 1000, 110, 2000), ncol = 2, byrow = TRUE)),
    list(ms_level = 1L, peaks = matrix(c(101, 1500), ncol = 2)),
    list(ms_level = 1L, peaks = matrix(c(102, 1200), ncol = 2))
  ), envir = .GlobalEnv)
  region <- asari_plot_mz_region(
    infile, tempfile(fileext = ".pdf"), c(0, 3), c(90, 120)
  )
  expect_true(file.exists(region))
})

test_that("PCA和相关图使用特征表样本列", {
  table <- .public_analysis_table(
    c("F1", "F2", "F3"), c(100, 200, 300), c(10, 20, 30),
    "S1", c(10, 20, 40)
  )
  table$S2 <- c(12, 19, 35)
  table$S3 <- c(8, 25, 50)
  path <- tempfile(fileext = ".tsv")
  utils::write.table(table, path, sep = "\t", row.names = FALSE, quote = FALSE)
  scores <- asari_pca(path, tempfile(fileext = ".pdf"))
  expect_equal(dim(scores), c(3L, 2L))
  correlation <- asari_plot_correlations(path, tempfile(fileext = ".pdf"))
  expect_true(file.exists(correlation))
})
