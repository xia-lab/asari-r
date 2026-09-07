# Each of the 9 defs in samples.py must have a clear R corresponding function.
test_that("all 9 Python samples defs have explicit R counterparts", {
  expected_functions <- c(
    "SimpleSample__init__",
    "SimpleSample_list_scan_numbers",
    "SimpleSample_get_mass_tracks_for_sample",
    "SimpleSample_get_masstracks_and_anchors",
    "SimpleSample_get_rt_calibration_records",
    "SimpleSample_extract_ms2",
    "SimpleSample__get_sample_data",
    "SimpleSample__retrieve_from_disk",
    "SimpleSample_load_intermediate"
  )

  # The number must be exactly 9, and each name resolves to a function.
  expect_length(expected_functions, 9L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# Create a minimal registry used by tests; fields are consistent with workflow.register_samples output.
make_sample_registry <- function(data_location = tempfile(fileext = ".json"),
                                 sample_data = list(
                                   list_mass_tracks = list(
                                     list(id_number = 0L, mz = 100)
                                   ),
                                   ms2_spectra = list()
                                 )) {
  list(
    input_file = "sample_1.mzML",
    name = "sample_1",
    sample_id = 0L,
    data_location = data_location,
    track_mzs = list(c(100, 0L)),
    max_scan_number = 3L,
    anchor_mz_pairs = list(c(0L, 1L), c(1L, 2L)),
    list_scan_numbers = 0:2,
    list_retention_time = c(0, 5, 10),
    sample_data = sample_data
  )
}

# Create a mutable experiment object to retain extract_ms2 updates to parameters.
make_sample_experiment <- function(output_dir = tempdir()) {
  experiment <- new.env(parent = emptyenv())
  experiment$parameters <- list(ms2_spectra_outdir = output_dir)
  experiment
}

# Verify that the constructor copies the registry field-by-field and initializes the Python object state correctly.
test_that("SimpleSample constructor preserves fields and mutable state", {
  registry <- make_sample_registry()
  experiment <- make_sample_experiment()
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory",
    mode = "neg",
    is_reference = TRUE
  )

  expect_s3_class(sample, "SimpleSample")
  expect_identical(sample$experiment, experiment)
  expect_equal(sample$mode, "neg")
  expect_true(sample$is_reference)
  expect_equal(sample$input_file, registry$input_file)
  expect_equal(sample$name, registry$name)
  expect_equal(sample$sample_id, 0L)
  expect_equal(sample$rt_numbers, 0:2)
  expect_equal(sample$list_scan_numbers, 0:2)
  expect_equal(sample$list_mass_tracks, registry$sample_data$list_mass_tracks)
  expect_equal(sample$`_mz_landmarks_`, c(0L, 1L, 2L))
  expect_equal(sample$rt_landmarks, list())
  expect_null(sample$rt_cal_dict)
  expect_null(sample$reverse_rt_cal_dict)
  expect_true(sample$is_rt_aligned)
  expect_null(sample$mz_calibration_function)
})

# Python properties are read from a private registry each time, rather than using a potentially out-of-date public copy.
test_that("list_scan_numbers property reads the current registry value", {
  registry <- list2env(make_sample_registry(), parent = emptyenv())
  sample <- SimpleSample__init__(registry, database_mode = "memory")

  registry$list_scan_numbers <- c(4L, 5L, 6L)

  expect_equal(SimpleSample_list_scan_numbers(sample), c(4L, 5L, 6L))
  expect_equal(sample$list_scan_numbers, c(4L, 5L, 6L))
  expect_error(
    sample$list_scan_numbers <- 7:9,
    "read-only property"
  )
})

# The memory mode should directly return the same mass track, and the static entry and instance entry results are consistent.
test_that("mass tracks are returned directly when already in memory", {
  registry <- make_sample_registry()
  sample <- SimpleSample__init__(registry, database_mode = "memory")

  expect_identical(
    SimpleSample_get_masstracks_and_anchors(sample),
    registry$sample_data$list_mass_tracks
  )
  expect_identical(
    SimpleSample_get_mass_tracks_for_sample(sample),
    registry$sample_data$list_mass_tracks
  )
  expect_identical(
    sample$get_masstracks_and_anchors(),
    registry$sample_data$list_mass_tracks
  )
})

# ondisk mode does not cache mass track, but loads JSON on demand through _get_sample_data.
test_that("mass tracks are loaded lazily in ondisk mode", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    '{"list_mass_tracks":[{"id_number":0,"mz":101.25}]}',
    json_path
  )
  registry <- make_sample_registry(data_location = json_path)
  sample <- SimpleSample__init__(registry, database_mode = "ondisk")

  expect_equal(sample$list_mass_tracks, list())
  tracks <- SimpleSample_get_masstracks_and_anchors(sample)
  expect_equal(tracks[[1L]]$id_number, 0L)
  expect_equal(tracks[[1L]]$mz, 101.25)
  expect_equal(SimpleSample__get_sample_data(sample),
               SimpleSample__retrieve_from_disk(sample))
})

# The calibration record must contain the four keys returned by the Python original function and reflect the latest state of the object.
test_that("RT calibration records contain current sample state", {
  sample <- SimpleSample__init__(make_sample_registry(), database_mode = "memory")
  sample$rt_landmarks <- c(10L, 20L)
  sample$reverse_rt_cal_dict <- c(`10` = 11L, `20` = 19L)

  records <- SimpleSample_get_rt_calibration_records(sample)

  expect_named(
    records,
    c("sample_id", "name", "rt_landmarks", "reverse_rt_cal_dict")
  )
  expect_equal(records$sample_id, 0L)
  expect_equal(records$name, "sample_1")
  expect_equal(records$rt_landmarks, c(10L, 20L))
  expect_equal(records$reverse_rt_cal_dict, c(`10` = 11L, `20` = 19L))
})

# Both normal JSON and the first JSON item in a ZIP should restore the same Python dict/list structure.
test_that("load_intermediate reads JSON and zipped JSON", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    '{"list_mass_tracks":[{"id_number":0,"mz":100.5}],"flag":true}',
    json_path
  )
  direct <- SimpleSample_load_intermediate(json_path)

  zip_dir <- tempfile("sample_zip_")
  dir.create(zip_dir)
  zip_json <- file.path(zip_dir, "sample.json")
  file.copy(json_path, zip_json)
  zip_path <- tempfile(fileext = ".zip")
  old_directory <- setwd(zip_dir)
  on.exit(setwd(old_directory), add = TRUE)
  suppressWarnings(utils::zip(zip_path, "sample.json"))
  zipped <- SimpleSample_load_intermediate(zip_path)

  expect_equal(direct, zipped)
  expect_equal(direct$list_mass_tracks[[1L]]$mz, 100.5)
  expect_true(direct$flag)
})

# The numpy tokens saved by json_tricks must be restored to R vectors rather than exposing the encoder internal dictionary.
test_that("load_intermediate decodes common json_tricks numpy markers", {
  skip_if_not_installed("jsonlite")
  json_path <- tempfile(fileext = ".json")
  writeLines(
    paste0(
      '{"intensity":{"__ndarray__":[1,2,3],',
      '"dtype":"int64","shape":[3]}}'
    ),
    json_path
  )

  decoded <- SimpleSample_load_intermediate(json_path)

  expect_equal(decoded$intensity, c(1L, 2L, 3L))
})

# Both Python pickle and pickle in ZIP recover common asari data types through conversion bridges.
test_that("load_intermediate reads Python pickle and zipped pickle", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle parity testing.")

  work_dir <- tempfile("pickle_test_")
  dir.create(work_dir)
  pickle_path <- file.path(work_dir, "sample.pickle")
  zip_path <- file.path(work_dir, "sample.zip")
  script_path <- file.path(work_dir, "make_pickle.py")
  writeLines(c(
    "import pickle, sys, zipfile",
    "path, zip_path = sys.argv[1:3]",
    "data = {'list_mass_tracks': [{'id_number': 0, 'mz': 100.25}], 'values': (1, 2)}",
    "with open(path, 'wb') as handle: pickle.dump(data, handle)",
    "with zipfile.ZipFile(zip_path, 'w') as archive: archive.write(path, 'sample.pickle')"
  ), script_path)
  output <- system2(
    python,
    shQuote(c(script_path, pickle_path, zip_path)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_null(attr(output, "status"))

  direct <- SimpleSample_load_intermediate(pickle_path)
  zipped <- SimpleSample_load_intermediate(zip_path)

  expect_equal(direct, zipped)
  expect_equal(direct$list_mass_tracks[[1L]]$mz, 100.25)
  expect_equal(unlist(direct$values), c(1L, 2L))
})

# Verify MS2 transformations, paths, and parameter side effects using the injectable Spectrum/save_spectra callback.
test_that("extract_ms2 converts spectra and calls the configured exporter", {
  captured <- new.env(parent = emptyenv())
  experiment <- make_sample_experiment(tempdir())
  experiment$Spectrum <- function(mz, intensities, metadata) {
    list(mz = mz, intensities = intensities, metadata = metadata)
  }
  experiment$save_spectra <- function(spectra, path, export_style) {
    captured$spectra <- spectra
    captured$path <- path
    captured$export_style <- export_style
    invisible(NULL)
  }
  ms2_spectrum <- list(
    scan_time_in_minutes = function() 1.5,
    mz = c(50, 75),
    intensity = c(1000, 500),
    precursor_mz = 150
  )
  registry <- make_sample_registry(sample_data = list(
    list_mass_tracks = list(),
    ms2_spectra = list(ms2_spectrum)
  ))
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory"
  )

  expect_null(SimpleSample_extract_ms2(sample, ".msp"))
  expect_equal(experiment$parameters$ms2_export_format, "msp")
  expect_equal(captured$path, file.path(tempdir(), "ms2_sample_1.msp"))
  expect_equal(captured$export_style, "matchms")
  expect_equal(captured$spectra[[1L]]$mz, c(50, 75))
  expect_equal(captured$spectra[[1L]]$intensities, c(1000, 500))
  expect_equal(captured$spectra[[1L]]$metadata$scan_time, 90)
  expect_equal(captured$spectra[[1L]]$metadata$origin, "sample_1")
  expect_equal(captured$spectra[[1L]]$metadata$precursor_mz, 150)
})

# Consistent with Python's try/except: missing MS2 export dependencies will only prompt, not throw an error.
test_that("extract_ms2 reports optional dependency failure without stopping", {
  experiment <- make_sample_experiment()
  registry <- make_sample_registry(sample_data = list(
    list_mass_tracks = list(),
    ms2_spectra = list()
  ))
  sample <- SimpleSample__init__(
    registry,
    experiment = experiment,
    database_mode = "memory"
  )

  expect_message(
    expect_null(SimpleSample_extract_ms2(sample)),
    "Error Extracting MS2 for: sample_1"
  )
  # An empty spectrum list will first update the format before failing at Python's save_spectra name resolution.
  expect_equal(experiment$parameters$ms2_export_format, "msp")
})

# Unknown extensions, unknown first ZIP items, and missing registry fields must all fail explicitly.
test_that("invalid intermediate files and incomplete registries fail clearly", {
  unknown_path <- tempfile(fileext = ".txt")
  writeLines("not an intermediate", unknown_path)
  expect_error(
    SimpleSample_load_intermediate(unknown_path),
    "Unknown file format"
  )
  expect_error(
    SimpleSample_load_intermediate(tempfile(fileext = ".json")),
    "Failed to load sample data"
  )
  expect_error(
    SimpleSample__init__(list()),
    "Missing registry field: input_file"
  )
})
