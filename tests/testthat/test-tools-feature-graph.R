source(file.path("R", "tools_feature_graph.R"), local = FALSE)

.feature_graph_table <- function() {
  file <- tempfile(fileext = ".tsv")
  data <- data.frame(
    id_number = c("F1", "F2", "F3"), mz = c(50, 60, 100),
    rtime = c(10, 10.2, 20), a = 1, b = 1, c = 1, d = 1, e = 1,
    f = 1, g = 1, h = 1, sample1 = c(10, 20, 0), sample2 = c(5, 10, 1),
    check.names = FALSE
  )
  utils::write.table(data, file, sep = "\t", row.names = FALSE, quote = FALSE)
  file
}

test_that("FeatureGraph建立完全图并按RT过滤", {
  graph <- FeatureGraph(.feature_graph_table())
  expect_s3_class(graph, "FeatureGraph")
  expect_equal(nrow(graph$graph$edges), 3L)
  expect_equal(FeatureGraph_metric(list(mz = 10, rtime = 2), list(mz = 12, rtime = 5)), c(2, 3))
  filtered <- graph$filter_graph(0.5)
  expect_equal(nrow(filtered$graph$edges), 1L)
  filtered$find_spectral_clusters()
  expect_length(filtered$clusters, 1L)
  expect_setequal(filtered$clusters[[1L]], c("F1", "F2"))
})

test_that("fragmentation spectrum选择峰数足够的最高强度样本", {
  graph <- FeatureGraph(.feature_graph_table())$filter_graph(0.5)
  spectra <- graph$extract_fragmentation_spectrum(TRUE, MIN_PEAKS_EXTRACTION = 2L)
  expect_length(spectra, 1L)
  expect_identical(spectra[[1L]]$metadata$sample, "sample1")
  expect_equal(max(spectra[[1L]]$intensities), 1)
})

test_that("注释传播到cluster内所有feature", {
  path <- .feature_graph_table()
  graph <- FeatureGraph(path)$filter_graph(0.5)
  graph$find_spectral_clusters()
  match <- list(list(
    extract = list(metadata = list(cluster_id = 0L)),
    library = list(metadata = list(compound_name = "C", inchikey = "K")),
    similarity = 0.9, match_peaks = 2L
  ))
  graph$map_annotations(match)
  expect_true(file.exists(sub(".tsv", "_annotated_gc_beta.tsv", path, fixed = TRUE)))
  expect_length(graph$df$annotations[[1L]], 3L)
})

test_that("FeatureGraph十个Python def均有R对应", {
  names <- c(
    "FeatureGraph__init__", "FeatureGraph_ftgraph_from_ft",
    "FeatureGraph_ftgraph_from_graph", "FeatureGraph_metric",
    "FeatureGraph_ft_to_graph", "FeatureGraph_graph_to_ft",
    "FeatureGraph_filter_graph", "FeatureGraph_find_spectral_clusters",
    "FeatureGraph_extract_fragmentation_spectrum", "FeatureGraph_map_annotations"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
