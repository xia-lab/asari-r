source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "asari_process.R"), local = FALSE)

test_that("公开入口只接受用户明确指定的mzML路径", {
  input_dir <- tempfile("asari-input-")
  dir.create(input_dir)
  file.create(file.path(input_dir, c("sample.mzML", "ignore.txt")))

  files <- .asari_collect_input_files(input_dir)
  expect_length(files, 1L)
  expect_identical(basename(files), "sample.mzML")
  expect_error(.asari_collect_input_files(tempfile()), "does not exist")
  expect_error(.asari_collect_input_files(file.path(input_dir, "ignore.txt")), "not a .mzML")
})

test_that("便捷参数覆盖asari默认配置并使用独立输出父目录", {
  output <- tempfile("asari-output-")
  prepared <- .asari_prepare_parameters(
    .asari_merge_parameters(list(min_peak_height = 5000)),
    output = output,
    project_name = "demo",
    mode = "neg",
    ppm = 8,
    multicores = 1,
    rt_align = FALSE,
    database_mode = "memory"
  )

  expect_identical(prepared$mode, "neg")
  expect_equal(prepared$mz_tolerance_ppm, 8)
  expect_identical(prepared$multicores, 1L)
  expect_false(prepared$rt_align_on)
  expect_identical(prepared$database_mode, "memory")
  expect_equal(prepared$min_peak_height, 5000)
  expect_true(dir.exists(output))
})

test_that("公开入口按Python规则补齐峰检测派生阈值", {
  prepared <- .asari_prepare_parameters(
    .asari_merge_parameters(list(min_peak_height = 300)),
    output = tempfile("asari-output-"),
    project_name = "demo",
    database_mode = "memory"
  )
  finalized <- .asari_finalize_peak_parameters(
    prepared,
    input_files = "explicit.mzML",
    override_names = "min_peak_height"
  )

  expect_identical(finalized$min_prominence_threshold, 99L)
  expect_equal(finalized$cal_min_peak_height, 3000)

  prepared$min_prominence_threshold <- 25
  overridden <- .asari_finalize_peak_parameters(
    prepared,
    input_files = "explicit.mzML",
    override_names = c("min_peak_height", "min_prominence_threshold")
  )
  expect_equal(overridden$min_prominence_threshold, 25)
})

test_that("参数错误会在读取数据前给出清楚提示", {
  expect_error(.asari_merge_parameters(c(a = 1)), "named list")
  expect_error(
    .asari_prepare_parameters(
      default_parameters(), tempfile(), "bad/name", database_mode = "memory"
    ),
    "without path separators"
  )
  expect_error(
    .asari_prepare_parameters(
      default_parameters(), tempfile(), "demo", ppm = 0, database_mode = "memory"
    ),
    "positive number"
  )
})

test_that("默认参数公开函数每次返回独立列表", {
  first <- asari_default_parameters()
  second <- asari_default_parameters()
  first$mode <- "neg"
  expect_identical(second$mode, "pos")
})
