# Locate the asari-r root directory from any working directory to avoid testthat not being able to find the R source code after switching directories.
arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
if (length(script_argument) == 1L) {
  script_file <- sub("^--file=", "", script_argument)
  package_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
} else if (file.exists("DESCRIPTION")) {
  package_root <- normalizePath(".", mustWork = TRUE)
} else {
  stop("Unable to locate the asari-r package root.", call. = FALSE)
}
setwd(package_root)

library(testthat)

# Direct the implicit plot device under test to a temporary file to avoid leaving Rplots.pdf in the source directory.
old_device <- getOption("device")
options(device = function(...) {
  grDevices::pdf(file = tempfile("asariR-test-", fileext = ".pdf"), ...)
})
on.exit(options(device = old_device), add = TRUE)

# First read all modules in the order in which the R package is loaded; some old tests only contain assertions and do not duplicate source target modules.
r_files <- sort(list.files("R", pattern = "[.]R$", full.names = TRUE))
invisible(lapply(r_files, source, local = .GlobalEnv))

# Manually register each test file so that the old file-by-file source testing and unified summary can work at the same time.
test_files <- sort(list.files("tests/testthat", pattern = "^test-.*[.]R$", full.names = TRUE))
reporter <- ListReporter$new()
with_reporter(reporter, {
  for (test_file in test_files) {
    reporter$start_file(basename(test_file))
    source(test_file, local = .GlobalEnv)
    reporter$end_file()
  }
})
results <- reporter$get_results()
results_frame <- as.data.frame(results)

# Outputs a stable total number, and returns a non-zero status on failure, which is convenient for Terminal and CI judgment.
passed <- sum(results_frame$passed)
failed <- sum(results_frame$failed) + sum(results_frame$error)
warnings <- sum(results_frame$warning)
skipped <- sum(results_frame$skipped)
cat(sprintf(
  "\nPASS %d | FAIL %d | WARN %d | SKIP %d\n",
  passed, failed, warnings, skipped
))
if (failed > 0L || warnings > 0L) {
  print(results_frame[
    results_frame$failed > 0L | results_frame$error |
      results_frame$warning > 0L,
    , drop = FALSE
  ])
  stop("The asariR test suite did not pass cleanly.", call. = FALSE)
}
