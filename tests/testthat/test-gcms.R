source(file.path("R", "mass_functions.R"), local = FALSE)
source(file.path("R", "tools_cosine.R"), local = FALSE)
source(file.path("R", "tools_match_features.R"), local = FALSE)
source(file.path("R", "tools_msp_parser.R"), local = FALSE)
source(file.path("R", "gcms.R"), local = FALSE)

.gcms_test_feature <- function(id, mz, ri, area, rtime = 60) {
  list(
    id = id, id_number = id, mz = mz, RI = ri, peak_area = area,
    rtime = rtime, cSelectivity = 1, goodness_fitting = 1,
    snr = 10, detection_counts = 2L, left_base = 1, right_base = 2,
    parent_masstrack_id = "M1"
  )
}

test_that("GC库重排、峰过滤和MSP输出对应Python", {
  records <- list(list(
    Name = "Compound", InChIKey = "KEY", RETENTIONTIME = "1200",
    ExactMass = "100", Formula = "C1", peaks = list(c(50, 10), c(60, 100))
  ))
  library <- reformat_gcms_lib(records, filter_factor = 20)
  expect_length(library, 1L)
  expect_s3_class(library[[1L]], "GC_lib_entry")
  expect_identical(library[[1L]]$base_peak, c(60, 100))
  expect_equal(nrow(library[[1L]]$peaks), 2L)
  expect_identical(filter_peaks_by_low_intensity_factor(records[[1L]]$peaks, 100, 20), records[[1L]]$peaks)

  outfile <- tempfile(fileext = ".msp")
  json_pseudospectra_to_msp(list(list(id = "P", RI = 1200, peaks = list(c(50, 10)))), outfile)
  expect_match(paste(readLines(outfile), collapse = "\n"), "Num Peaks: 1", fixed = TRUE)
})

test_that("seed伪谱和RI函数保留严格窗口规则", {
  features <- list(
    .gcms_test_feature("F1", 50, 1000, 100),
    .gcms_test_feature("F2", 60, 1020, 20)
  )
  frame <- data.frame(S1 = c(10, 5), S2 = c(20, 10), row.names = c("F1", "F2"))
  pseudo <- get_seeded_pseudospectrum(
    "seed", 1000, features[[1L]], features, frame,
    max_ri_delta = 100, low_peak_filter_factor = 1000,
    feature_distance_filter = NULL
  )
  expect_s3_class(pseudo, "PseudoSpectrum")
  expect_equal(pseudo$num_features, 2L)
  expect_equal(ri_penalty_function(0), 1)
  expect_equal(ri_penalty_function(101), 0)
  expect_length(find_entries_in_rtwindow(1000, list(GC_lib_entry("c", "", "", 1030, 0, "", 50, matrix(c(50, 1), ncol = 2), c(50, 1), "")), 30), 0L)
})

test_that("库评分整理和数组序列化形成一致结构", {
  features <- list(
    .gcms_test_feature("F1", 50, 1000, 100),
    .gcms_test_feature("F2", 60, 1001, 50)
  )
  dict <- stats::setNames(features, c("F1", "F2"))
  frame <- data.frame(S1 = c(10, 5), S2 = c(20, 10), row.names = c("F1", "F2"))
  lib <- GC_lib_entry(
    "C1", "KEY", "Compound", 1000, 100, "C1",
    c(50, 60), matrix(c(50, 100, 60, 50), ncol = 2, byrow = TRUE),
    c(50, 100), ""
  )
  scored <- batch_lib_search_score(list(lib), features, dict, frame, corr_cutoff = 0)
  expect_length(scored, 1L)
  curated <- curate_batch_lib_search_result(scored, score_cutoff_cosine = 0, score_cutoff_entropy = 0)
  expect_length(curated[[1L]], 1L)
  serialized <- serialize_annotated_empCpds(curated[[1L]])
  expect_true(is.list(serialized[[1L]]$peaks_in_lib))
})

test_that("重复注释清理合并余弦和熵分数", {
  base <- list(
    id = "E", inchikey = "K", quant_ion = "F1", features = c("F1"),
    score = 0.8, peaks_as_features = matrix(c(50, 1), ncol = 2),
    peaks_in_lib = matrix(c(50, 1), ncol = 2)
  )
  feature <- list(feature = "F1", inchikey = "K", quant_ion = "F1", score = 0.8)
  cleaned <- cleanup_anno_empcpds_features(list(
    list_empCpds_cosine = list(base), feature_anno_list_cosine = list(feature),
    list_empCpds_entropy = list(base), feature_anno_list_entropy = list(feature),
    core_features = "F1"
  ))
  expect_equal(cleaned[[1L]]$K[[1L]]$score_cosine, 0.8)
  expect_equal(cleaned[[1L]]$K[[1L]]$score_entropy, 0.8)
  expect_equal(cleaned[[2L]]$F1[[1L]]$score_entropy, 0.8)
})

test_that("gcms全部45个Python def和两个NamedTuple类均有对应", {
  functions <- c(
    "load_gcms_dbfile", "reformat_gcms_lib", "json_pseudospectra_to_msp",
    "designate_base_peak", "filter_peaks_by_intensity_factor",
    "filter_peaks_by_low_intensity_factor", "filter_features_by_low_intensity_factor",
    "filter_against_libentry", "find_entries_in_rtwindow",
    "read_fit_KovatsIndex_rtime", "append_kovats_index",
    "port_pseudospectrum_to_json", "serialize_annotated_empCpds",
    "ri_penalty_function", "filter_peaks_by_penalized_distance",
    "find_all_matches_centurion_indexed_list", "get_matched_features_per_cpd",
    "distill_correlated_features", "get_seeded_pseudospectrum",
    "get_spaced_top_features", "iterative_build_pseudospectra_by_penalizeddistance",
    "have_basepeak_molecularion", "batch_lib_search_score",
    "curate_batch_lib_search_result", "write_tsv_feature_anno",
    "write_tsv_empCpd_anno", "get_clusters_by_ri_hcl", "extend_cluster",
    "format_fcluster_to_pseudospectrum", "iterative_build_pseudospectra_by_hcl",
    "batch_lib_search_by_basepeaks", "export_feature_annotation_bybasepeaksearch",
    "group_pseudospectra_from_features", "reverse_spec_searches",
    "export_feature_annotations", "export_feature_annotation_details",
    "iterative_reverse_annotation", "is_same_match", "is_same_matched_feature",
    "merge_epd_with_list", "merge_feature_with_list", "format_arrays_to_lists",
    "cleanup_anno_empcpds_features", "append_anno_json", "export_jsonanno_tsv"
  )
  expect_length(functions, 45L)
  expect_true(all(vapply(functions, function(name) is.function(get(name)), FALSE)))
  expect_true(is.function(PseudoSpectrum))
  expect_true(is.function(GC_lib_entry))
})
