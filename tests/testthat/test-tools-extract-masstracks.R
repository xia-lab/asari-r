source(file.path("R", "tools_extract_masstracks.R"), local = FALSE)

test_that("read_project_dir_masstracks按Python子串规则筛选", {
  directory <- tempfile()
  dir.create(directory)
  file.create(file.path(directory, c("a.mzML", "b.mzML.gz", "c.txt")))
  expect_setequal(
    basename(read_project_dir_masstracks(directory)),
    c("a.mzML", "b.mzML.gz")
  )
})

test_that("process_single_file写出轨迹统计", {
  old <- getOption("asariR.simple_sample_from_file")
  on.exit(options(asariR.simple_sample_from_file = old), add = TRUE)
  sample <- new.env(parent = emptyenv())
  sample$get_mass_tracks_ <- function(...) {
    sample$list_mass_tracks <- list(list(id_number = 0L, mz = 100, intensity = c(0, 2, 5)))
  }
  options(asariR.simple_sample_from_file = function(path) sample)
  outdir <- tempfile()
  dir.create(outdir)
  process_single_file("sample.mzML", outdir)
  lines <- readLines(file.path(outdir, "sample_massTrack.tsv"))
  expect_identical(lines[[1L]], "id_number\tmz\tnumber_scans\tmax_intensity\tmedian_intensity")
  expect_identical(lines[[2L]], "0\t100\t3\t5\t3")
})

test_that("两个Python def均有R对应", {
  expect_true(is.function(read_project_dir_masstracks))
  expect_true(is.function(process_single_file))
})
