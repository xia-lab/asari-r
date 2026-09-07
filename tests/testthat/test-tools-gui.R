source(file.path("R", "tools_gui.R"), local = FALSE)

test_that("TextRedirector写入widget并支持flush", {
  widget <- new.env(parent = emptyenv())
  widget$text <- ""
  widget$insert <- function(where, value) widget$text <- paste0(widget$text, value)
  widget$see <- function(where) invisible(NULL)
  redirector <- TextRedirector(widget)
  redirector$write("hello")
  redirector$flush()
  expect_identical(widget$text, "hello")
})

test_that("参数编辑和选项回调保留类型", {
  old_edit <- getOption("asariR.edit_parameters")
  old_select <- getOption("asariR.select_option")
  on.exit(options(asariR.edit_parameters = old_edit, asariR.select_option = old_select), add = TRUE)
  options(asariR.edit_parameters = function(values, types) { values$n <- "3"; values })
  edited <- create_ui(list(flag = TRUE, n = 1L, x = NULL, unsupported = list(1)))
  expect_true(edited$flag)
  expect_equal(edited$n, 3L)
  expect_null(edited$x)
  options(asariR.select_option = function(options, key) options[[2L]])
  expect_identical(create_selection_ui(c("a", "b"), list(), "run")$run, "b")
})

test_that("gui全部15个Python def含内部回调均有对应", {
  names <- c(
    "run_program", "select_directory", "TextRedirector__init__",
    "TextRedirector_write", "TextRedirector_flush", "run_program_thread",
    "start_program_gui", "show_disclaimer", ".gui_on_accept", ".gui_on_decline",
    "create_ui", ".gui_run_callback", "create_selection_ui", ".gui_set_value", "main_gui"
  )
  expect_length(names, 15L)
  expect_true(all(vapply(names, function(name) is.function(get(name)), FALSE)))
})
