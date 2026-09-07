source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "main.R"), local = FALSE)

test_that("main的12个Python def都有明确R对应", {
  names <- c(
    "process", "analyze", "annotate", "join", "viz",
    "update_peak_detection_params", "update_params_from_CLI", "__debug_print",
    "initialize_parameters", "build_parser", "run_asari", "main"
  )
  expect_length(names, 12L)
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})

test_that("build_parser解析长短选项和数值类型", {
  args <- build_parser(c(
    "process", "-m", "neg", "--ppm=8", "-i", ".",
    "--autoheight", "True", "--workflow", "LC", "-c", "1"
  ))
  expect_identical(args$run, "process")
  expect_identical(args$mode, "neg")
  expect_identical(args$ppm, 8L)
  expect_identical(args$multicores, 1L)
  expect_identical(args$autoheight, "True")
})

test_that("CLI覆盖参数并保留Python的ppm字段行为", {
  args <- build_parser(c(
    "process", "--mode", "neg", "--ppm", "7", "--output", tempdir(),
    "--autoheight", "False", "--anno", "True", "--peak_area", "auc"
  ))
  parameters <- update_params_from_CLI(default_parameters(), args)
  expect_identical(parameters$mode, "neg")
  expect_identical(parameters$ppm, 7L)
  expect_identical(parameters$mz_tolerance_ppm, 5)
  expect_false(parameters$autoheight)
  expect_true(parameters$anno)
  expect_identical(parameters$peak_area, "auc")
  expect_identical(parameters$run, "process")
})

test_that("峰检测派生阈值与Python相同", {
  parameters <- default_parameters()
  parameters$autoheight <- FALSE
  updated <- update_peak_detection_params(
    parameters,
    list(min_peak_height = 300, min_prominence_threshold = NULL,
         cal_min_peak_height = NULL, min_intensity_threshold = 25)
  )
  expect_equal(updated$min_peak_height, 300)
  expect_identical(updated$min_prominence_threshold, 99L)
  expect_equal(updated$cal_min_peak_height, 3000)
  expect_equal(updated$min_intensity_threshold, 25)
})

test_that("初始化写入版本和14位时间戳", {
  initialized <- initialize_parameters(default_parameters(), list())
  expect_identical(initialized$asari_version, "1.17.0")
  expect_match(initialized$timestamp, "^[0-9]{8}-[0-9]{6}$")
})

test_that("list_workflows现在返回可直接使用的工作流", {
  output <- capture.output(result <- run_asari(list(run = "list_workflows"), list()))
  expect_true(is.data.frame(result))
  expect_true(any(grepl("asari_process", output, fixed = TRUE)))
})
