source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "tools_cosine.R"), local = FALSE)
source(file.path("R", "tools_msp_parser.R"), local = FALSE)
source(file.path("R", "tools_feature_graph.R"), local = FALSE)
source(file.path("R", "tools_gc_annotation.R"), local = FALSE)

test_that("EI_MS_Library读取用户MSP并标准化谱图", {
  file <- tempfile(fileext = ".msp")
  writeLines(c("Name: Test", "Num Peaks: 2", "50 10", "60 20"), file)
  library <- EI_MS_Library(file, multicores = 1L)
  expect_s3_class(library, "EI_MS_Library")
  expect_length(library$library, 1L)
  expect_equal(max(library$library[[1L]]$intensities), 1)
  expect_identical(library$library_meta$Parser, "matchms")
})

test_that("wrapped_cosine保留job和后端评分", {
  old <- getOption("asariR.cosine_greedy_pair")
  on.exit(options(asariR.cosine_greedy_pair = old), add = TRUE)
  options(asariR.cosine_greedy_pair = function(x, y) c(0.8, 3L))
  job <- list(list(id = "a"), list(id = "b"))
  result <- wrapped_cosine(job)
  expect_identical(result[[1L]], job)
  expect_equal(result[[2L]], c(0.8, 3))
})

test_that("GC_Annotation按模型写回registry RI", {
  annotation <- GC_Annotation(
    parameters = list(),
    sample_registry = list(S = list(list_retention_time = c(10, 20)))
  )
  annotation$RI_map <- list(existing = TRUE)
  annotation$RI_models <- list(R = function(x) x * 100)
  annotation$convert_to_RI(list(S = "R"))
  expect_equal(annotation$sample_registry$S$list_retention_index, c(1000, 2000))
})

test_that("两个类的十个Python def均有R对应", {
  names <- c(
    "GC_Annotation_populate_RI_lookup", "GC_Annotation_convert_to_RI",
    "GC_Annotation_annotate_GC", "EI_MS_Library__init__",
    "EI_MS_Library_load_library", "EI_MS_Library_retrieve_library_meta",
    "EI_MS_Library_load_library_manifest", "EI_MS_Library_annotate_gc_feature_table",
    "EI_MS_Library_annotate_gc_feature_table_with_library", "wrapped_cosine"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
