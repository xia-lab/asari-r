source(file.path("R", "default_parameters.R"), local = FALSE)
source(file.path("R", "utils.R"), local = FALSE)
source(file.path("R", "main.R"), local = FALSE)
source(file.path("R", "command_line.R"), local = FALSE)

test_that("command_line与main使用同一入口", {
  output <- capture.output(result <- command_line_main(c("list_workflows")))
  expect_identical(result$run, "list_workflows")
  expect_true(any(grepl("Hello from Asari", output, fixed = TRUE)))
})
