# Each of the 10 defs in utils.py must have a clear R corresponding function.
test_that("all 10 Python utils defs have explicit R counterparts", {
  expected_functions <- c(
    "NpEncoder_default",
    "bulk_process",
    "download_and_unzip_to_pkg_resources",
    "download_and_unzip",
    "validate_mzml_file",
    "build_boolean_dict",
    "sizeof_fmt",
    "checksum_file",
    "wait_with_pbar",
    "get_ionization_mode_mzml"
  )

  expect_length(expected_functions, 10L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# NpEncoder preserves class identity and converts scalar and ndarray corresponding structures to base R types.
test_that("NpEncoder default converts numeric and array values", {
  encoder <- NpEncoder()

  expect_s3_class(encoder, "NpEncoder")
  expect_identical(encoder$default(3L), 3L)
  expect_equal(encoder$default(3.5), 3.5)
  expect_equal(
    encoder$default(matrix(1:6, nrow = 2L, byrow = TRUE)),
    list(as.list(1:3), as.list(4:6))
  )
  expect_equal(encoder$default(c(1, 2)), list(1, 2))
  expect_error(encoder$default(list(a = 1)), "not JSON serializable")
})

# bulk_process must maintain input order and pass each argument to command as the only argument.
test_that("bulk_process applies command in argument order", {
  result <- bulk_process(
    function(value) value * value,
    as.list(1:5),
    jobs_per_worker = 1
  )

  expect_equal(unlist(result), c(1, 4, 9, 16, 25))
  expect_equal(bulk_process(identity, list(), jobs_per_worker = 1), list())
  expect_error(
    bulk_process("not a function", list(1), jobs_per_worker = 1),
    "command must be a function"
  )
})

# Create a local ZIP to avoid downloading and testing access to the external network.
make_utils_zip <- function() {
  source_dir <- tempfile("utils_zip_source_")
  dir.create(source_dir)
  writeLines("payload", file.path(source_dir, "payload.txt"))
  archive <- tempfile(fileext = ".zip")
  old_directory <- setwd(source_dir)
  on.exit(setwd(old_directory), add = TRUE)
  suppressWarnings(utils::zip(archive, "payload.txt"))
  archive
}

# download_and_unzip should download and restore the file contents from the ZIP.
test_that("download_and_unzip extracts a local archive", {
  archive <- make_utils_zip()
  target <- tempfile("utils_extract_")

  expect_message(
    expect_null(download_and_unzip(paste0("file://", archive), target)),
    "Extracted to"
  )
  expect_true(file.exists(file.path(target, "payload.txt")))
  expect_equal(readLines(file.path(target, "payload.txt")), "payload")
})

# The package resource version is decompressed to the subdir under the parent of the package directory according to Python logic.
test_that("download_and_unzip_to_pkg_resources resolves package parent", {
  archive <- make_utils_zip()
  package_parent <- tempfile("fake_library_")
  package_path <- file.path(package_parent, "fakePackage")
  dir.create(package_path, recursive = TRUE)
  subdir <- "downloaded_data"
  expected_path <- file.path(package_parent, subdir, "payload.txt")

  expect_message(
    expect_null(download_and_unzip_to_pkg_resources(
      paste0("file://", archive), package_path, subdir
    )),
    "Extracted to"
  )
  expect_true(file.exists(expected_path))
})

# Invalid or non-existent mzML must return FALSE rather than propagating a parsing exception to the caller.
test_that("validate_mzml_file rejects invalid files", {
  invalid <- tempfile(fileext = ".mzML")
  writeLines("not mzML", invalid)

  expect_false(validate_mzml_file(invalid))
  expect_false(validate_mzml_file(tempfile(fileext = ".mzML")))
})

# All ten Boolean input forms supported by Python must be mapped to correct logical values.
test_that("build_boolean_dict reproduces all Python keys", {
  booleans <- build_boolean_dict()

  expect_named(
    booleans,
    c("T", "F", "1", "0", "True", "False", "TRUE", "FALSE", "true", "false")
  )
  expect_equal(
    unlist(booleans, use.names = FALSE),
    c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
})

# sizeof_fmt's unit rounding, negative numbers, and custom suffixes are consistent with Python format.
test_that("sizeof_fmt uses Python IEC units and formatting", {
  expect_equal(sizeof_fmt(0), "0.0B")
  expect_equal(sizeof_fmt(1023), "1023.0B")
  expect_equal(sizeof_fmt(1024), "1.0KiB")
  expect_equal(sizeof_fmt(-1024^2), "-1.0MiB")
  expect_equal(sizeof_fmt(1536, suffix = "bytes"), "1.5Kibytes")
  expect_equal(sizeof_fmt(1024^8), "1.0YiB")
})

# checksum_file returns standard lowercase MD5 and rejects directories or non-existent paths.
test_that("checksum_file computes MD5", {
  path <- tempfile()
  writeChar("hello", path, eos = NULL)

  expect_equal(checksum_file(path), "5d41402abc4b2a76b9719d911017c592")
  expect_error(checksum_file(tempdir()), "regular file")
  expect_error(checksum_file(tempfile()), "regular file")
  expect_error(checksum_file(path, chunksize = 0), "Invalid chunksize")
})

# wait=0 completes immediately as does Python range(0), illegal counts should fail explicitly.
test_that("wait_with_pbar handles zero and invalid waits", {
  expect_null(wait_with_pbar(0))
  expect_error(wait_with_pbar(-1), "non-negative integer")
  expect_error(wait_with_pbar(1.5), "non-negative integer")
})

# The polarity header covers positive, negative, mixed, limit and unset modes.
test_that("get_ionization_mode_mzml identifies scan polarity", {
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(1, 1)))),
    "pos"
  )
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(0, 0)))),
    "neg"
  )
  expect_equal(
    get_ionization_mode_mzml(list(header = data.frame(polarity = c(1, 0)))),
    "mixed"
  )
  expect_equal(
    get_ionization_mode_mzml(
      list(header = data.frame(polarity = c(1, 0))), limit = 0
    ),
    "pos"
  )
  expect_error(
    get_ionization_mode_mzml(list(header = data.frame(polarity = -1))),
    "neither positive nor negative"
  )
})

# Keep the ppm conversion compatibility functions already used by the original R project.
test_that("ppm_to_mz_tolerance converts ppm to absolute m/z tolerance", {
  expect_equal(ppm_to_mz_tolerance(100, 5), 0.0005)
  expect_equal(ppm_to_mz_tolerance(800, 5), 0.004)
})
