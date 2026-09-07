test_that("mass_functions全15个Python def包含两个内部函数都有对应", {
  functions <- c(
    "flatten_tuplelist", "check_close_mzs", "calculate_selectivity", "__sel__",
    "mass_paired_mapping", "complete_mass_paired_mapping", "all_mass_paired_mapping",
    "_find_all_mzmatches_centurion_indexed_list",
    "mass_paired_mapping_with_correction", "landmark_guided_mapping",
    "mass_functions_bin_by_median", "gap_divide_mz_cluster", "__divide_by_largest_gap__",
    "identify_mass_peaks", "nn_cluster_by_mz_seeds"
  )
  expect_length(functions, 15L)
  expect_true(all(vapply(functions, function(name) is.function(get(name)), FALSE)))
})

test_that("flatten_tuplelist flattens tuple pairs and removes duplicates", {
  tuples <- list(c(1, 2), c(2, 3), c(3, 4))

  expect_equal(flatten_tuplelist(tuples), c(1, 2, 3, 4))
})

test_that("flatten_tuplelist validates its tuple input", {
  expect_null(flatten_tuplelist(list()))
  expect_error(flatten_tuplelist(c(1, 2)), "must be a list")
  expect_error(flatten_tuplelist(list(1, c(2, 3))), "at least two")
})

test_that("check_close_mzs returns one-based adjacent close pairs", {
  mzs <- c(100.0000, 100.0003, 100.0020, 100.0022)

  expect_equal(
    check_close_mzs(mzs, mz_tolerance_ppm = 5),
    list(c(2L, 1L), c(4L, 3L))
  )
})

test_that("check_close_mzs checks adjacent values rather than all pairs", {
  mzs <- c(100.0000, 100.0002, 100.0003)

  expect_equal(
    check_close_mzs(mzs, mz_tolerance_ppm = 5),
    list(c(2L, 1L), c(3L, 2L))
  )
})

test_that("check_close_mzs handles short input and validates sorting", {
  expect_equal(check_close_mzs(numeric()), list())
  expect_equal(check_close_mzs(100), list())
  expect_error(check_close_mzs(c(100.1, 100)), "ascending order")
})

test_that("check_close_mzs uses a strict tolerance comparison", {
  upper_mz <- 100
  exact_lower_boundary <- upper_mz * (1 - 5e-6)

  expect_equal(
    check_close_mzs(c(exact_lower_boundary, upper_mz), 5),
    list()
  )
  expect_equal(
    check_close_mzs(c(exact_lower_boundary + 1e-7, upper_mz), 5),
    list(c(2L, 1L))
  )
})

test_that("calculate_selectivity matches the Python four-value calculation", {
  mzs <- c(100.0000, 100.0005, 100.0010, 100.0015)
  expected <- c(
    0.54657191622213974,
    0.34549779899513522,
    0.34549632904743122,
    0.54656787982816335
  )

  expect_equal(calculate_selectivity(mzs, std_ppm = 5), expected, tolerance = 1e-14)
})

test_that("calculate_selectivity handles interior values and cutoff branches", {
  expect_equal(
    calculate_selectivity(c(100, 100.000001, 100.001, 100.02, 100.2)),
    c(0, 0, 0.74741079425419588, 1, 1),
    tolerance = 1e-14
  )
})

test_that("calculate_selectivity保留Python对非升序输入的实际行为", {
  expect_error(calculate_selectivity(c(100, 101, 102)), "more than three")
  expect_equal(
    calculate_selectivity(c(100, 100.2, 100.1, 100.3)),
    c(1, 0, 0, 1)
  )
  expect_error(calculate_selectivity(100:103, std_ppm = 0), "positive")
})

test_that("mass_paired_mapping matches unsorted lists and preserves positions", {
  result <- mass_paired_mapping(
    c(200.0000, 100.0000),
    c(200.0008, 100.0003),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(2L, 2L), c(1L, 1L)))
  expect_equal(
    result$ratio_deltas,
    c(2.999990999984481e-6, 3.999984000054678e-6),
    tolerance = 1e-16
  )
})

test_that("mass_paired_mapping leaves an ambiguous candidate unmatched", {
  result <- mass_paired_mapping(100, c(100.0002, 100.0003), std_ppm = 5)

  expect_equal(result$mapped, list())
  expect_equal(result$ratio_deltas, numeric())
})

test_that("mass_paired_mapping keeps the list2 minus list1 shift sign", {
  result <- mass_paired_mapping(100.0003, 100, std_ppm = 5)

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(
    result$ratio_deltas,
    -2.999990999984481e-6,
    tolerance = 1e-16
  )
})

test_that("mass_paired_mapping handles empty lists and validates inputs", {
  expect_equal(
    mass_paired_mapping(numeric(), numeric()),
    list(mapped = list(), ratio_deltas = numeric())
  )
  expect_error(mass_paired_mapping("100", 100), "numeric vectors")
  expect_error(mass_paired_mapping(100, 100, std_ppm = -1), "non-negative")
})

test_that("complete_mass_paired_mapping chooses the closest candidate", {
  result <- complete_mass_paired_mapping(
    c(100.0000, 150.0000, 200.0000),
    c(100.0002, 100.0003, 150.0004, 300.0000),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 3L)))
  expect_equal(result$list1_unmapped, 3L)
  expect_equal(result$list2_unmapped, c(2L, 4L))
})

test_that("complete_mass_paired_mapping resolves a shared-list2 conflict", {
  result <- complete_mass_paired_mapping(
    c(99.9998, 100.0001),
    100,
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(2L, 1L)))
  expect_equal(result$list1_unmapped, 1L)
  expect_equal(result$list2_unmapped, integer())
})

test_that("complete_mass_paired_mapping returns all unmatched positions", {
  expect_equal(
    complete_mass_paired_mapping(c(100, 200), c(300, 400)),
    list(
      mapped = list(),
      list1_unmapped = c(1L, 2L),
      list2_unmapped = c(1L, 2L)
    )
  )
  expect_equal(
    complete_mass_paired_mapping(numeric(), 100),
    list(
      mapped = list(),
      list1_unmapped = integer(),
      list2_unmapped = 1L
    )
  )
})

test_that("complete_mass_paired_mapping validates inputs", {
  expect_error(
    complete_mass_paired_mapping("100", 100),
    "numeric vectors"
  )
  expect_error(
    complete_mass_paired_mapping(100, 100, std_ppm = -1),
    "non-negative"
  )
})

test_that("all_mass_paired_mapping keeps every candidate within tolerance", {
  result <- all_mass_paired_mapping(
    c(100.0000, 100.0004, 200.0000),
    c(100.0002, 300.0000),
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 1L)))
  expect_equal(result$list1_unmapped, 3L)
  expect_equal(result$list2_unmapped, 2L)
})

test_that("all_mass_paired_mapping searches adjacent Centurion buckets", {
  result <- all_mass_paired_mapping(99.9999, 100.0001, std_ppm = 5)

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(result$list1_unmapped, integer())
  expect_equal(result$list2_unmapped, integer())
})

test_that("all_mass_paired_mapping uses a strict ppm tolerance", {
  result <- all_mass_paired_mapping(
    c(100.00049, 100.00051),
    100,
    std_ppm = 5
  )

  expect_equal(result$mapped, list(c(1L, 1L)))
  expect_equal(result$list1_unmapped, 2L)
  expect_equal(result$list2_unmapped, integer())
})

test_that("all_mass_paired_mapping returns unmatched positions", {
  expect_equal(
    all_mass_paired_mapping(c(100, 200), c(300, 400)),
    list(
      mapped = list(),
      list1_unmapped = c(1L, 2L),
      list2_unmapped = c(1L, 2L)
    )
  )
  expect_equal(
    all_mass_paired_mapping(numeric(), numeric()),
    list(
      mapped = list(),
      list1_unmapped = integer(),
      list2_unmapped = integer()
    )
  )
})

test_that("all_mass_paired_mapping validates inputs", {
  expect_error(all_mass_paired_mapping("100", 100), "numeric vectors")
  expect_error(
    all_mass_paired_mapping(100, 100, std_ppm = -1),
    "non-negative"
  )
})

test_that("mass mapping correction can recover a new positive-shift match", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200, 300),
    c(100.0004, 200.0008, 300.0018),
    std_ppm = 5,
    correction_tolerance_ppm = 1
  )

  expect_equal(
    result$mapped,
    list(c(1L, 1L), c(2L, 2L), c(3L, 3L))
  )
  expect_equal(
    result$correction_ratio,
    3.999984000054678e-6,
    tolerance = 1e-16
  )
})

test_that("mass mapping correction leaves a sub-threshold shift unchanged", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200),
    c(100.00005, 200.0001),
    correction_tolerance_ppm = 1
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 2L)))
  expect_equal(
    result$correction_ratio,
    4.999997500167233e-7,
    tolerance = 1e-17
  )
})

test_that("mass mapping correction preserves Python negative-shift behavior", {
  result <- mass_paired_mapping_with_correction(
    c(100, 200),
    c(99.9997, 199.9994),
    correction_tolerance_ppm = 1
  )

  expect_equal(result$mapped, list(c(1L, 1L), c(2L, 2L)))
  expect_equal(
    result$correction_ratio,
    -2.999999999957481e-6,
    tolerance = 1e-16
  )
})

test_that("mass mapping correction returns NaN when there are no anchors", {
  result <- mass_paired_mapping_with_correction(100, 200)

  expect_equal(result$mapped, list())
  expect_true(is.nan(result$correction_ratio))
})

test_that("mass mapping correction validates inputs", {
  expect_error(
    mass_paired_mapping_with_correction("100", 100),
    "numeric vectors"
  )
  expect_error(
    mass_paired_mapping_with_correction(100, 100, std_ppm = -1),
    "std_ppm"
  )
  expect_error(
    mass_paired_mapping_with_correction(
      100,
      100,
      correction_tolerance_ppm = -1
    ),
    "correction_tolerance_ppm"
  )
})

test_that("landmark guided mapping corrects and extends the reference", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200, 300),
      c(1L, 2L),
      c(100.0004, 200.0008, 300.0018, 400.0016),
      c(1L, 2L, 4L),
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(
    result$new_reference_mzlist,
    c(
      100.00000000079999,
      200.00000000159997,
      300.0003000012,
      400.00000000639994
    ),
    tolerance = 1e-12
  )
  expect_equal(result$new_reference_map2, c(1L, 2L, 3L, 4L))
  expect_equal(result$REF_landmarks, c(1L, 2L, 4L))
  expect_equal(
    result$correction_ratio,
    3.999984000054678e-6,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping uses NA for an absent sample track", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200),
      1L,
      c(100.00005, 300),
      1L,
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(result$new_reference_mzlist, c(100.000025, 200, 300))
  expect_equal(result$new_reference_map2, c(1L, NA_integer_, 2L))
  expect_equal(result$REF_landmarks, 1L)
  expect_equal(
    result$correction_ratio,
    4.999997500167233e-7,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping corrects a negative mass shift", {
  result <- suppressMessages(
    landmark_guided_mapping(
      c(100, 200, 300),
      c(1L, 2L),
      c(99.9997, 199.9994, 299.9991),
      c(1L, 2L),
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_equal(result$new_reference_map2, c(1L, 2L, 3L))
  expect_equal(
    result$correction_ratio,
    -2.999999999957481e-6,
    tolerance = 1e-16
  )
})

test_that("landmark guided mapping leaves correction unset with too few anchors", {
  result <- suppressMessages(
    landmark_guided_mapping(
      100 * seq_len(10),
      1L,
      c(100.0003, 1100),
      1L,
      std_ppm = 5,
      correction_tolerance_ppm = 1
    )
  )

  expect_null(result$correction_ratio)
  expect_equal(result$new_reference_map2[[1L]], 1L)
  expect_equal(tail(result$new_reference_map2, 1L), 2L)
})

test_that("landmark guided mapping validates landmark positions", {
  expect_error(
    landmark_guided_mapping(100, 0, 100, 1L),
    "REF_landmarks"
  )
  expect_error(
    landmark_guided_mapping(100, 1L, 100, 2L),
    "SM_landmarks"
  )
})

test_that("bin_by_median follows the moving median boundary", {
  tuples <- list(
    list(100.0, "a"),
    list(100.2, "b"),
    list(100.9, "c"),
    list(102.0, "d")
  )

  expect_equal(
    bin_by_median(tuples, function(value) 0.5),
    list(list("a", "b"), list("c"), list("d"))
  )
})

test_that("bin_by_median uses a strict tolerance comparison", {
  tuples <- list(list(100, "a"), list(100.5, "b"))

  expect_equal(
    bin_by_median(tuples, function(value) 0.5),
    list(list("a"), list("b"))
  )
})

test_that("gap_divide_mz_cluster splits at the first largest gap", {
  tuples <- list(
    c(100.0, 1, 10),
    c(100.25, 2, 20),
    c(101.25, 3, 30),
    c(102.25, 4, 40),
    c(102.5, 5, 50)
  )

  expect_equal(
    gap_divide_mz_cluster(tuples, mz_tolerance = 0.1),
    list(tuples[1:2], tuples[3:5])
  )
})

test_that("nearest integer filter matches SciPy even and odd windows", {
  expect_equal(
    .uniform_filter1d_nearest_integer(1:7, 2),
    c(1L, 1L, 2L, 3L, 4L, 5L, 6L)
  )
  expect_equal(
    .uniform_filter1d_nearest_integer(1:7, 5),
    c(1L, 2L, 3L, 4L, 5L, 5L, 6L)
  )
})

test_that("peak selection matches SciPy plateau and distance behavior", {
  expect_equal(.find_peaks_with_distance(c(0, 1, 1, 0), 1), 2L)
  expect_equal(.find_peaks_with_distance(c(0, 1, 1, 1, 0), 1), 3L)
  expect_equal(
    .find_peaks_with_distance(c(0, 2, 0, 2, 0), 3),
    4L
  )
  expect_equal(
    .find_peaks_with_distance(c(0, 1, 0, 1, 0, 1, 0), 3),
    c(2L, 6L)
  )
})

test_that("identify_mass_peaks reproduces SciPy integer smoothing", {
  mzs <- c(
    rep(100.0000, 5),
    rep(100.0001, 3),
    rep(100.0015, 6),
    rep(100.0016, 2)
  )
  tuples <- lapply(seq_along(mzs), function(ii) c(mzs[[ii]], ii, 1))

  expect_equal(
    identify_mass_peaks(tuples, mz_tolerance = 0.0005),
    100.0014
  )
})

test_that("identify_mass_peaks sorts only when presorted is false", {
  mzs <- rep(c(100.0015, 100.0000, 100.0016, 100.0001), 4)
  tuples <- lapply(seq_along(mzs), function(ii) c(mzs[[ii]], ii, 1))

  expect_equal(
    identify_mass_peaks(
      tuples,
      mz_tolerance = 0.0005,
      presorted = FALSE
    ),
    numeric()
  )
})

test_that("nn_cluster_by_mz_seeds assigns each tuple to the nearest seed", {
  mzs <- c(
    rep(100.000, 3),
    rep(100.001, 8),
    rep(100.002, 10),
    rep(100.003, 4)
  )
  tuples <- lapply(seq_along(mzs), function(ii) c(mzs[[ii]], ii, 1))
  clusters <- nn_cluster_by_mz_seeds(
    tuples,
    mz_tolerance = 0.0005,
    presorted = TRUE
  )

  expect_length(clusters, 2L)
  expect_equal(lengths(clusters), c(11L, 14L))
  expect_equal(
    vapply(clusters[[1L]], `[[`, numeric(1), 1L),
    mzs[1:11]
  )
  expect_equal(
    vapply(clusters[[2L]], `[[`, numeric(1), 1L),
    mzs[12:25]
  )
})

test_that("nn_cluster_by_mz_seeds falls back to the largest gap", {
  tuples <- list(c(100.0000, 1, 1), c(100.0020, 2, 1))

  expect_equal(
    nn_cluster_by_mz_seeds(
      tuples,
      mz_tolerance = 0.0005,
      presorted = TRUE
    ),
    list(tuples[1], tuples[2])
  )
})

test_that("remaining mass functions validate required inputs", {
  expect_error(bin_by_median(list(), function(value) 1), "non-empty")
  expect_error(gap_divide_mz_cluster(list(c(100, 1)), 1), "at least two")
  expect_error(identify_mass_peaks(list(), 0.0005), "non-empty")
  expect_error(
    identify_mass_peaks(list(c(100, 1)), 0.00001),
    "at least 0.0001"
  )
})
