source(file.path("R", "gcms.R"), local = FALSE)
source(file.path("R", "dashboard.R"), local = FALSE)

.dashboard_fixture <- function() {
  peak <- list(
    id_number = "F1", parent_masstrack_id = 0L, snr = 10,
    goodness_fitting = 0.95, cSelectivity = 0.96, height = 100,
    peak_area = 1000, mz = 100, rtime_left_base = 1,
    rtime = 2, rtime_right_base = 3
  )
  cmap <- list(
    `_number_of_samples_` = 1L, rt_length = 3L,
    dict_scan_rtime = list(`0` = 0, `1` = 1, `2` = 2),
    list_mass_tracks = list(list(id_number = 0L, mz = 100, intensity = c(0, 10, 0))),
    rt_records = list(list(name = "S1", reverse_rt_cal_dict = list(`0` = 0L, `2` = 2L)))
  )
  list(peak = peak, cmap = cmap)
}

test_that("epd_convert和track_to_peaks建立双向索引", {
  fixture <- .dashboard_fixture()
  converted <- epd_convert(list(E1 = list(MS1_pseudo_Spectra = list(fixture$peak))))
  expect_named(converted[[1L]], "F1")
  expect_identical(converted[[2L]]$E1$MS1_pseudo_Spectra, "F1")
  expect_identical(track_to_peaks(converted[[1L]])$`0`, "F1")
})

test_that("轨迹搜索和RT对齐数据对应Python", {
  fixture <- .dashboard_fixture()
  expect_equal(find_track_by_mz(fixture$cmap, 0:2, 100.05), 0L)
  alignment <- prepare_rt_alignment(fixture$cmap)
  expect_equal(alignment$S1, c(0, 0, 0))
  plot <- cmapplot_mass_tracks(fixture$cmap, 0:2, "blue", 0L)
  expect_equal(plot$data$intensity, c(0, 10, 0))
})

test_that("dashboard返回可调用浏览状态", {
  fixture <- .dashboard_fixture()
  epd <- list(E1 = list(MS1_pseudo_Spectra = list(fixture$peak)))
  table <- data.frame(mz = 100, rtime = 2, peak_area = 1000, snr = 10, goodness_fitting = 0.95, cSelectivity = 0.96)
  state <- dashboard(list(outdir = "x", project_name = "p"), fixture$cmap, epd, table)
  expect_equal(state$default_feature_number, 1L)
  expect_match(state$feature_info_by_feature_id(1), "id_number", fixed = TRUE)
  expect_equal(nrow(state$features_by_rt_with_lines(2, 0.5)), 1L)
})

test_that("dashboard全部18个Python def含六个内部函数均有对应", {
  names <- c(
    "epd_convert", "read_project", "plot_xic", "cmapplot_mass_tracks",
    "convert_dict_html", "convert_dict_markdown", "track_to_peaks",
    "find_track_by_mz", "find_a_good_peak", "prepare_rt_alignment",
    "get_summary_panel", "dashboard", ".dashboard_feature_info_by_feature_id",
    ".dashboard_cmapplot_track_by_feature_id", ".dashboard_get_features_by_mz",
    ".dashboard_track_info_by_mz", ".dashboard_cmapplot_track_by_mz",
    ".dashboard_features_by_rt_with_lines"
  )
  expect_length(names, 18L)
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
