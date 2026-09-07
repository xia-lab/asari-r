source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "tools_mzml_converter.R"), local = FALSE)

test_that("converter构造命令并批量替换输入输出", {
  directory <- tempfile(); dir.create(directory)
  file.create(file.path(directory, "ThermoRawFileParser.exe"))
  .mzMLconverter_installed$`1.4.5` <- directory
  old_probe <- getOption("asariR.converter_probe")
  old_run <- getOption("asariR.converter_run")
  on.exit(options(asariR.converter_probe = old_probe, asariR.converter_run = old_run), add = TRUE)
  options(asariR.converter_probe = function(command) length(command) == 1L)
  seen <- list()
  options(asariR.converter_run = function(command) { seen[[length(seen) + 1L]] <<- command; 0L })
  converter <- mzMLconverter(multicores = 1L)
  converter$bulk_convert(file.path(directory, "sample.raw"))
  expect_match(paste(seen[[1L]], collapse = " "), "sample.raw", fixed = TRUE)
  expect_match(paste(seen[[1L]], collapse = " "), directory, fixed = TRUE)
})

test_that("mzml_converter六个Python def含内部函数均有对应", {
  names <- c(
    "mzMLconverter__init__", "mzMLconverter__conversion_command_template",
    "mzMLconverter__determine_executable", "mzMLconverter_install_converter",
    "mzMLconverter_uninstall_converter", "mzMLconverter_bulk_convert"
  )
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
