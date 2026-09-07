# Each of the 12 defs in workflow.py must have a clear R corresponding function.
test_that("all 12 Python workflow defs have explicit R counterparts", {
  expected_functions <- c(
    "workflow_setup",
    "process_project",
    "read_project_dir",
    "read_project_file",
    "register_samples",
    "create_export_folders",
    "remove_intermediate_pickles",
    "make_iter_parameters",
    "batch_EIC_from_samples_",
    "single_sample_EICs_",
    "process_xics",
    "get_mz_list"
  )

  expect_length(expected_functions, 12L)
  expect_true(all(vapply(expected_functions, function(name) {
    exists(name, mode = "function", inherits = TRUE)
  }, logical(1))))
})

# Generate complete parameters required for workflow testing; use independent output prefixes each time to avoid timestamp collisions.
make_workflow_parameters <- function(database_mode = "memory",
                                     workflow = "DIMS") {
  list(
    outdir = tempfile("workflow_project_"),
    project_name = "asari",
    time_stamp_for_dir = "fixed",
    database_mode = database_mode,
    project_sample_number_small = 10,
    workflow = workflow,
    anno = FALSE,
    mode = "pos",
    mz_tolerance_ppm = 5,
    check_isotope_ratio = FALSE,
    min_intensity_threshold = 100,
    min_timepoints = 3,
    min_peak_height = 1000,
    multicores = 1,
    compress = FALSE,
    reuse_intermediates = NULL,
    rt_align_on = FALSE,
    drop_unaligned_samples = FALSE,
    output_feature_table = "Feature_table.tsv",
    mass_grid_mapping = "_mass_grid_mapping.csv",
    correction_tolerance_ppm = 1,
    cal_min_peak_height = 1000,
    peak_number_rt_calibration = 2,
    num_lowess_iterations = 3,
    max_retention_shift = NULL,
    debug_rtime_align = FALSE,
    peak_area = "sum",
    min_peak_ratio = 0.001,
    min_intensity = 100,
    wlen = 5,
    signal_noise_ratio = 2,
    gaussian_shape = 0.5,
    min_prominence_threshold = 100
  )
}

# Return deterministic five mass track so that the real MassGrid and selectivity paths meet the input conditions.
make_workflow_extractor <- function(acquisition_time = 1700000000) {
  force(acquisition_time)
  function(infile, mz_tolerance_ppm, min_intensity,
           min_timepoints, min_peak_height) {
    offset <- if (grepl("sample_b", basename(infile), fixed = TRUE)) 0.0001 else 0
    mzs <- 100:104 + offset
    intensity <- c(0L, 10L, 100L, 10L, 0L)
    list(
      rt_numbers = 0:4,
      rt_times = seq(0, 20, 5),
      tracks = lapply(seq_along(mzs), function(ii) {
        list(mzs[[ii]], as.integer(intensity / ii))
      }),
      ms2_spectra = list(),
      acquisition_time = acquisition_time
    )
  }
}

# The anchor lookup mock returns two valid Python 0-based track pairs.
workflow_anchor_finder <- function(list_mass_tracks, mz_tolerance_ppm) {
  stopifnot(length(list_mass_tracks) == 5L, mz_tolerance_ppm == 5)
  list(c(0L, 1L), c(2L, 3L))
}

# Sample registration must use consecutive 0-based IDs and preserve the order of incoming files.
test_that("register_samples creates Python style registry", {
  files <- c("a.mzML", "b.mzML")
  registry <- register_samples(files)

  expect_named(registry, c("0", "1"))
  expect_equal(registry[["0"]]$sample_id, 0L)
  expect_equal(registry[["1"]]$sample_id, 1L)
  expect_equal(
    unname(vapply(registry, `[[`, character(1), "input_file")),
    files
  )
})

# Directory reading is filtered by file name substring; after the project list is read, it is converted to absolute path and blanks are removed.
test_that("project readers filter mzML paths", {
  directory <- tempfile("mzml_dir_")
  dir.create(directory)
  file.create(file.path(directory, c("a.mzML", "b.txt", "c.mzML.gz")))

  from_dir <- read_project_dir(directory)
  expect_equal(
    sort(basename(from_dir)),
    c("a.mzML", "c.mzML.gz")
  )

  project_file <- tempfile(fileext = ".txt")
  writeLines(c(" a.mzML ", "ignore.txt", "sub/b.mzML"), project_file)
  from_file <- read_project_file(project_file)
  expect_equal(
    basename(from_file),
    c("a.mzML", "b.mzML")
  )
  expect_true(all(startsWith(from_file, normalizePath(getwd()))))
})

# Update outdir/export/tmp_pickle_dir when creating the project directory and refuse to overwrite targets with the same name.
test_that("create_export_folders creates Python directory layout", {
  parameters <- make_workflow_parameters()
  updated <- create_export_folders(parameters, "12345")

  expect_true(dir.exists(updated$outdir))
  expect_true(dir.exists(updated$export_outdir))
  expect_true(dir.exists(updated$tmp_pickle_dir))
  expect_equal(updated$outdir, paste(
    parameters$outdir, parameters$project_name, "12345", sep = "_"
  ))

  already <- updated
  expect_message(
    returned <- create_export_folders(already, "other"),
    "already exist"
  )
  expect_equal(returned$outdir, updated$outdir)
})

# A temporary pickle cleanup only deletes the specified project directory and refuses deletion in reuse mode.
test_that("remove_intermediate_pickles removes only temporary pickle files", {
  parameters <- create_export_folders(make_workflow_parameters(), "cleanup")
  files <- file.path(parameters$tmp_pickle_dir, c("a.pickle", "b.pickle"))
  file.create(files)

  expect_null(remove_intermediate_pickles(parameters))
  expect_false(dir.exists(parameters$tmp_pickle_dir))
  expect_true(dir.exists(parameters$outdir))

  parameters$reuse_intermediates <- "old_project"
  expect_error(
    remove_intermediate_pickles(parameters),
    "Cannot remove"
  )
})

# make_iter_parameters should generate a four-tuple of sample_id, input, output, and the same parameter object.
test_that("make_iter_parameters builds one extraction job per sample", {
  parameters <- create_export_folders(make_workflow_parameters(), "jobs")
  registry <- register_samples(c("a.mzML", "b.mzML"))
  jobs <- make_iter_parameters(registry, parameters)

  expect_length(jobs, 2L)
  expect_equal(jobs[[1L]][[1L]], 0L)
  expect_equal(jobs[[1L]][[2L]], "a.mzML")
  expect_equal(jobs[[1L]][[3L]], file.path(
    parameters$outdir, "pickle", "a.pickle"
  ))
  expect_identical(jobs[[1L]][[4L]], parameters)
})

# Memory mode single sample extraction should return a Python 12-item tuple and full sample_data.
test_that("single sample extraction returns memory data and anchor metadata", {
  parameters <- make_workflow_parameters(database_mode = "memory")
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder
  outfile <- tempfile(fileext = ".pickle")
  result <- single_sample_EICs_(list(0L, "sample_a.mzML", outfile, parameters))
  tuple <- result[["0"]]

  expect_length(tuple, 12L)
  expect_equal(tuple[[1L]], "passed")
  expect_equal(tuple[[2L]], "passed")
  expect_equal(tuple[[4L]], 4L)
  expect_equal(tuple[[5L]], 0:4)
  expect_equal(tuple[[6L]], seq(0, 20, 5))
  expect_length(tuple[[7L]], 5L)
  expect_equal(tuple[[8L]], 2L)
  expect_equal(tuple[[9L]], list(c(0L, 1L), c(2L, 3L)))
  expect_equal(tuple[[10L]], 1700000000)
  expect_length(tuple[[11L]]$list_mass_tracks, 5L)
  expect_false(file.exists(outfile))
})

# ondisk mode must write real pickles that can be read by both samples.R and Python.
test_that("single sample extraction writes readable ondisk pickle", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle testing.")
  parameters <- make_workflow_parameters(database_mode = "ondisk")
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder
  outfile <- tempfile(fileext = ".pickle")

  tuple <- single_sample_EICs_(
    list(3L, "sample_a.mzML", outfile, parameters)
  )[["3"]]
  restored <- SimpleSample_load_intermediate(outfile)

  expect_equal(tuple[[1L]], "passed")
  expect_true(file.exists(outfile))
  expect_equal(restored$sample_id, 3L)
  expect_length(restored$list_mass_tracks, 5L)
  expect_equal(
    lapply(restored$anchor_mz_pairs, unlist, use.names = FALSE),
    list(c(0L, 1L), c(2L, 3L))
  )
  expect_equal(tuple[[11L]], list())
})

# Any fetch exceptions are converted to a 12-item tuple with status failed, without interrupting other samples.
test_that("single sample extraction converts errors to failed status", {
  parameters <- make_workflow_parameters()
  parameters$extract_mass_tracks <- function(...) stop("broken mzML")
  result <- NULL

  expect_message(
    result <- single_sample_EICs_(
      list(2L, "broken.mzML", "unused.pickle", parameters)
    ),
    "Failed to extract broken.mzML"
  )
  tuple <- result[["2"]]
  expect_length(tuple, 12L)
  expect_equal(tuple[[1L]], "failed")
  expect_equal(tuple[[2L]], "failed")
  expect_null(tuple[[3L]])
  expect_equal(tuple[[12L]], FALSE)
})

# Batch extraction merges each single sample result by sample ID without losing the order and quantity.
test_that("batch extraction merges all sample results", {
  parameters <- create_export_folders(make_workflow_parameters(), "batch")
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder
  registry <- register_samples(c("sample_a.mzML", "sample_b.mzML"))

  result <- batch_EIC_from_samples_(registry, parameters)

  expect_named(result, c("0", "1"))
  expect_true(all(vapply(result, function(tuple) {
    identical(tuple[[1L]], "passed")
  }, logical(1))))
})

# workflow_setup corrects the misalignment of Python tuples and constructs an experiment that can actually run DIMS.
test_that("workflow setup maps extraction tuple fields and runs real DIMS", {
  parameters <- make_workflow_parameters(database_mode = "auto")
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder
  experiment <- workflow_setup(
    c("sample_a.mzML", "sample_b.mzML"),
    parameters
  )

  expect_s3_class(experiment, "ext_Experiment")
  expect_equal(experiment$parameters$database_mode, "memory")
  expect_equal(experiment$sample_registry[["0"]]$acquisition_time, 1700000000)
  expect_true(is.list(experiment$sample_registry[["0"]]$sample_data))
  expect_equal(experiment$sample_registry[["0"]]$sparsified, FALSE)
  expect_equal(experiment$sample_registry[["0"]]$name, "sample_a")

  expect_null(experiment$process_all_DIMS())
  expect_equal(nrow(experiment$CMAP$MassGrid), 5L)
  expect_equal(nrow(experiment$CMAP$FeatureTable), 5L)
  expect_length(experiment$all_samples, 2L)
})

# Unsupported workflows will still complete setup and print Python Error 100, without calling the processing method by mistake.
test_that("process_project reports unsupported workflow without crashing", {
  parameters <- make_workflow_parameters(
    database_mode = "memory", workflow = "UNKNOWN"
  )
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder

  expect_message(
    expect_null(process_project(c("sample_a.mzML"), parameters)),
    "Workflow not supported. Error 100"
  )
})

# process_xics fixes Python's double register_samples bug and prints pickle.
test_that("process_xics stores one pickle per input mzML", {
  skip_if_not_installed("jsonlite")
  python <- tryCatch(.samples_find_python(), error = function(error) "")
  skip_if(!nzchar(python), "Python is required for pickle testing.")
  parameters <- make_workflow_parameters(database_mode = "memory")
  base_outdir <- parameters$outdir
  parameters$extract_mass_tracks <- make_workflow_extractor()
  parameters$find_mzdiff_pairs <- workflow_anchor_finder

  expect_message(
    expect_null(process_xics(c("sample_a.mzML", "sample_b.mzML"), parameters)),
    "XICs were stored"
  )
  projects <- Sys.glob(paste0(base_outdir, "_asari_*"))
  expect_length(projects, 1L)
  pickles <- list.files(file.path(projects, "pickle"), pattern = "\\.pickle$")
  expect_equal(sort(pickles), c("sample_a.pickle", "sample_b.pickle"))
})

# get_mz_list supports tab or comma delimitation in the first column, and rejects content that cannot be converted into numeric values.
test_that("get_mz_list reads target mz values from first column", {
  target_file <- tempfile(fileext = ".tsv")
  writeLines(c(
    "mz\tname",
    "100.125\tA",
    "200.25,comment\tB"
  ), target_file)
  expect_equal(get_mz_list(target_file), c(100.125, 200.25))

  invalid_file <- tempfile(fileext = ".tsv")
  writeLines(c("mz", "not-a-number"), invalid_file)
  expect_error(get_mz_list(invalid_file), "Invalid m/z value")
})
