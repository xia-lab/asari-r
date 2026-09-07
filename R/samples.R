# Corresponds to Python asari/samples.py: lightweight sample objects and intermediate data reading.
#
# Data structure convention:
# 1. Python’s mutable SimpleSample object uses environment in R to ensure that field modifications are retained.
# 2. Python dict corresponds to a named list or environment, None corresponds to NULL, and an empty list corresponds to list().
# 3. Business numbers such as scan and track maintain the original Python values; this module does not perform 0/1-based conversion without authorization.
# 4. Python pickle is not a native R format, so it is converted to JSON and read through the configurable Python interpreter.

# Reads a field from a named list or environment, returning the given default value if missing.
.samples_get <- function(object, name, default = NULL) {
  if (is.environment(object)) {
    if (exists(name, envir = object, inherits = FALSE)) {
      return(get(name, envir = object, inherits = FALSE))
    }
    return(default)
  }
  if (is.list(object) && name %in% names(object)) return(object[[name]])
  default
}

# Read required fields according to Python registry['field'] semantics; if missing, an error will be reported immediately.
.samples_registry_get_required <- function(registry, name) {
  sentinel <- new.env(parent = emptyenv())
  value <- .samples_get(registry, name, sentinel)
  if (identical(value, sentinel)) {
    stop("Missing registry field: ", name, call. = FALSE)
  }
  value
}

# Determine whether the mass track container is equivalent to a non-empty list in Python.
.samples_has_items <- function(value) {
  !is.null(value) && length(value) > 0L
}

# Modify the experiment's parameters field; the experiment object must be mutable to preserve Python's side effects.
.samples_set_experiment_parameter <- function(experiment, name, value) {
  if (!is.environment(experiment)) {
    stop("experiment must be an environment to update parameters.", call. = FALSE)
  }
  parameters <- .samples_get(experiment, "parameters", list())
  parameters[[name]] <- value
  experiment$parameters <- parameters
  invisible(value)
}

# Recursively restore the tag structures used by json_tricks for numpy ndarray and set.
.samples_decode_json_tricks <- function(value) {
  if (!is.list(value)) return(value)

  if ("__ndarray__" %in% names(value)) {
    data <- unlist(
      lapply(value$`__ndarray__`, .samples_decode_json_tricks),
      use.names = FALSE
    )
    dtype <- .samples_get(value, "dtype", "")
    if (grepl("^(u?int)", dtype)) data <- as.integer(data)
    if (grepl("^bool", dtype)) data <- as.logical(data)

    # Numpy saves multidimensional arrays row-wise; dimension padding is first reversed and then transposed back to Python axis order.
    shape <- as.integer(unlist(.samples_get(value, "shape", length(data))))
    if (length(shape) <= 1L) return(data)
    reversed <- array(data, dim = rev(shape))
    return(aperm(reversed, rev(seq_along(shape))))
  }
  if ("__set__" %in% names(value)) {
    return(lapply(value$`__set__`, .samples_decode_json_tricks))
  }

  lapply(value, .samples_decode_json_tricks)
}

# Parse JSON text into an R list and restore common special types for Python json_tricks.
.samples_parse_json <- function(text) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "Reading JSON intermediate files requires the optional jsonlite package.",
      call. = FALSE
    )
  }
  parsed <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  .samples_decode_json_tricks(parsed)
}

# Check whether the file has a ZIP PK file header to avoid warnings when calling unzip on ordinary files.
.samples_is_zipfile <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- readBin(connection, what = "raw", n = 4L)
  length(signature) >= 2L &&
    identical(as.integer(signature[1:2]), c(0x50L, 0x4bL))
}

# Read the raw bytes of the first item of the ZIP; the order is consistent with Python z.namelist()[0].
.samples_read_first_zip_entry <- function(path) {
  entries <- utils::unzip(path, list = TRUE)
  if (nrow(entries) == 0L) {
    return(list(name = NULL, raw = raw()))
  }

  entry_name <- entries$Name[[1L]]
  entry_size <- as.integer(entries$Length[[1L]])
  connection <- unz(path, entry_name, open = "rb")
  on.exit(close(connection), add = TRUE)
  list(
    name = entry_name,
    raw = readBin(connection, what = "raw", n = entry_size)
  )
}

# Locates the Python used by the pickle transformation; the user can specify this explicitly via option or environment variables.
.samples_find_python <- function() {
  configured <- getOption("asariR.python", "")
  if (nzchar(configured)) return(configured)

  configured <- Sys.getenv("ASARIR_PYTHON", unset = "")
  if (nzchar(configured)) return(configured)

  candidates <- c(
    file.path(".venv", "bin", "python"),
    file.path("..", ".venv", "bin", "python"),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Reading Python pickle files requires a Python interpreter. ",
      "Set options(asariR.python = '/path/to/python').",
      call. = FALSE
    )
  }
  existing[[1L]]
}

# Decode pickle through the Python standard library and convert it into JSON that R can read.
.samples_read_pickle <- function(data_location, zip_entry = NULL) {
  python <- .samples_find_python()
  script_path <- tempfile(fileext = ".py")
  output_path <- tempfile(fileext = ".json")
  on.exit(unlink(c(script_path, output_path)), add = TRUE)

  # Converters handle common numpy arrays, numpy scalars, and spectral objects with __dict__ common in asari.
  script <- c(
    "import json, pickle, sys, zipfile",
    "source, entry, output = sys.argv[1:4]",
    "if entry:",
    "    with zipfile.ZipFile(source, 'r') as archive:",
    "        obj = pickle.loads(archive.read(entry))",
    "else:",
    "    with open(source, 'rb') as handle:",
    "        obj = pickle.load(handle)",
    "def convert(value):",
    "    try:",
    "        import numpy as np",
    "        if isinstance(value, np.ndarray): return value.tolist()",
    "        if isinstance(value, np.generic): return value.item()",
    "    except ImportError:",
    "        pass",
    "    if isinstance(value, (set, tuple)): return list(value)",
    "    if isinstance(value, bytes): return {'__bytes__': list(value)}",
    "    if hasattr(value, '__dict__'):",
    "        payload = dict(value.__dict__)",
    "        scan_time = getattr(value, 'scan_time_in_minutes', None)",
    "        if callable(scan_time): payload['scan_time_in_minutes'] = scan_time()",
    "        for name in ('mz', 'intensity', 'precursor_mz'):",
    "            try: payload[name] = getattr(value, name)",
    "            except Exception: pass",
    "        return payload",
    "    raise TypeError(f'Unsupported pickle value: {type(value)!r}')",
    "with open(output, 'w', encoding='utf-8') as handle:",
    "    json.dump(obj, handle, default=convert, allow_nan=True)"
  )
  writeLines(script, script_path, useBytes = TRUE)

  status <- system2(
    python,
    # system2 leaves it to the shell to parse arguments, so file paths containing spaces and empty strings are explicitly quoted.
    args = shQuote(c(
      script_path,
      data_location,
      if (is.null(zip_entry)) "" else zip_entry,
      output_path
    )),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop(
      "Failed to decode Python pickle: ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }
  .samples_parse_json(paste(readLines(output_path, warn = FALSE), collapse = "\n"))
}

# Find optional callbacks that the Python version of extract_ms2 depends on, either from the experiment object or the global environment.
.samples_resolve_callback <- function(self, name) {
  callback <- .samples_get(self$experiment, name)
  if (is.function(callback)) return(callback)
  get0(name, envir = .GlobalEnv, mode = "function", inherits = TRUE)
}

# Corresponds to SimpleSample.__init__: Creates a variable sample object according to the registry.
SimpleSample__init__ <- function(registry = list(),
                                 experiment = NULL,
                                 database_mode = "ondisk",
                                 mode = "pos",
                                 is_reference = FALSE) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("SimpleSample", "environment")

  # Save construction parameters and private registry; field names correspond to Python objects one by one.
  self$experiment <- experiment
  self$mode <- mode
  self$database_mode <- database_mode
  self$is_reference <- is_reference
  self$`.__registry` <- registry

  self$input_file <- .samples_registry_get_required(registry, "input_file")
  self$name <- .samples_registry_get_required(registry, "name")
  self$sample_id <- .samples_registry_get_required(registry, "sample_id")
  self$data_location <- .samples_registry_get_required(registry, "data_location")
  self$track_mzs <- .samples_registry_get_required(registry, "track_mzs")
  self$max_scan_number <- .samples_registry_get_required(
    registry, "max_scan_number"
  )
  self$anchor_mz_pairs <- .samples_registry_get_required(
    registry, "anchor_mz_pairs"
  )
  self$rt_numbers <- .samples_registry_get_required(
    registry, "list_scan_numbers"
  )
  self$list_retention_time <- .samples_registry_get_required(
    registry, "list_retention_time"
  )

  # Memory mode holds mass track directly; ondisk mode only reads when needed.
  if (identical(database_mode, "memory")) {
    sample_data <- .samples_registry_get_required(registry, "sample_data")
    self$list_mass_tracks <- .samples_registry_get_required(
      sample_data, "list_mass_tracks"
    )
  } else {
    self$list_mass_tracks <- list()
  }

  # The m/z landmark flatten logic multiplexes mass_functions.R and initializes the RT calibration state.
  self$`_mz_landmarks_` <- flatten_tuplelist(self$anchor_mz_pairs)
  self$rt_landmarks <- list()
  self$rt_cal_dict <- NULL
  self$reverse_rt_cal_dict <- NULL
  self$is_rt_aligned <- is_reference
  self$mz_calibration_function <- NULL

  # Active binding corresponds to Python @property: the private registry is re-read on each access.
  makeActiveBinding(
    "list_scan_numbers",
    local({
      sample <- self
      function(value) {
        if (!missing(value)) {
          stop("list_scan_numbers is a read-only property.", call. = FALSE)
        }
        SimpleSample_list_scan_numbers(sample)
      }
    }),
    self
  )

  # Install method closures for environment so that constructors.R can be called as a Python object.
  self$get_mass_tracks_for_sample <- function(sample) {
    SimpleSample_get_mass_tracks_for_sample(sample)
  }
  self$get_masstracks_and_anchors <- function() {
    SimpleSample_get_masstracks_and_anchors(self)
  }
  self$get_rt_calibration_records <- function() {
    SimpleSample_get_rt_calibration_records(self)
  }
  self$extract_ms2 <- function(export_format = "msp") {
    SimpleSample_extract_ms2(self, export_format)
  }
  self$`_get_sample_data` <- function() SimpleSample__get_sample_data(self)
  self$`_retrieve_from_disk` <- function() SimpleSample__retrieve_from_disk(self)

  self
}

# Corresponds to @property SimpleSample.list_scan_numbers: always read from the original registry.
SimpleSample_list_scan_numbers <- function(self) {
  .samples_registry_get_required(self$`.__registry`, "list_scan_numbers")
}

# Corresponds to @staticmethod SimpleSample.get_mass_tracks_for_sample.
SimpleSample_get_mass_tracks_for_sample <- function(sample) {
  SimpleSample_get_masstracks_and_anchors(sample)
}

# Corresponds to SimpleSample.get_masstracks_and_anchors: priority is given to using memory, otherwise disk reading is required.
SimpleSample_get_masstracks_and_anchors <- function(self) {
  if (.samples_has_items(self$list_mass_tracks)) {
    return(self$list_mass_tracks)
  }
  sample_data <- SimpleSample__get_sample_data(self)
  .samples_registry_get_required(sample_data, "list_mass_tracks")
}

# Corresponds to SimpleSample.get_rt_calibration_records: exports the four fields required for RT calibration.
SimpleSample_get_rt_calibration_records <- function(self) {
  list(
    sample_id = self$sample_id,
    name = self$name,
    rt_landmarks = self$rt_landmarks,
    reverse_rt_cal_dict = self$reverse_rt_cal_dict
  )
}

# Corresponds to SimpleSample.extract_ms2: converts the MS2 spectrum and calls the optional saver.
SimpleSample_extract_ms2 <- function(self, export_format = "msp") {
  tryCatch({
    # The Python version reads from the registry in memory mode, and reads from intermediate files in other modes.
    if (identical(self$database_mode, "memory")) {
      sample_data <- .samples_registry_get_required(
        self$`.__registry`, "sample_data"
      )
    } else {
      sample_data <- SimpleSample__get_sample_data(self)
    }
    ms2_data <- .samples_registry_get_required(sample_data, "ms2_spectra")

    Spectrum <- .samples_resolve_callback(self, "Spectrum")
    save_spectra <- .samples_resolve_callback(self, "save_spectra")

    spectra <- lapply(ms2_data, function(spec) {
      # Python only parses Spectrum when actually traversing the spectrum; this constructor is not required for an empty spectrum list.
      if (!is.function(Spectrum)) stop("Spectrum callback is unavailable.")
      scan_time <- .samples_get(spec, "scan_time_in_minutes")
      if (is.null(scan_time)) stop("MS2 spectrum has no scan time.")
      if (is.function(scan_time)) scan_time <- scan_time()
      rtime <- scan_time * 60
      precursor_mz <- tryCatch(
        .samples_get(spec, "precursor_mz", NULL),
        error = function(error) NULL
      )

      # Like Python zip(spec.mz, spec.intensity), only the common length between the two is retained.
      mzs <- .samples_get(spec, "mz", numeric())
      intensities <- .samples_get(spec, "intensity", numeric())
      common_length <- min(length(mzs), length(intensities))
      if (common_length == 0L) {
        mzs <- intensities <- numeric()
      } else {
        mzs <- as.numeric(mzs[seq_len(common_length)])
        intensities <- as.numeric(intensities[seq_len(common_length)])
      }

      Spectrum(
        mz = mzs,
        intensities = intensities,
        metadata = list(
          scan_time = rtime,
          origin = self$name,
          precursor_mz = precursor_mz
        )
      )
    })

    # Python removes only the leading dot and writes the final format back to the experimental parameters.
    if (!is.character(export_format) || length(export_format) != 1L ||
        !nzchar(export_format)) {
      stop("export_format must contain at least one character.")
    }
    if (startsWith(export_format, ".")) {
      export_format <- substring(export_format, 2L)
    }
    .samples_set_experiment_parameter(
      self$experiment, "ms2_export_format", export_format
    )
    output_dir <- .samples_registry_get_required(
      self$experiment$parameters, "ms2_spectra_outdir"
    )
    path <- file.path(
      output_dir,
      sprintf("ms2_%s.%s", self$name, export_format)
    )
    if (!is.function(save_spectra)) {
      stop("save_spectra callback is unavailable.")
    }
    save_spectra(spectra, path, export_style = "matchms")
    invisible(NULL)
  }, error = function(error) {
    # Python catches all exceptions and only reports the current sample without blocking the main process.
    message("Error Extracting MS2 for: ", self$name)
    invisible(NULL)
  })
}

# Corresponds to SimpleSample._get_sample_data: Preserves the wrapping layer required for future database schema extensions.
SimpleSample__get_sample_data <- function(self) {
  SimpleSample_load_intermediate(self$data_location)
}

# Corresponds to SimpleSample._retrieve_from_disk: reads the intermediate file from the data location of the current sample.
SimpleSample__retrieve_from_disk <- function(self) {
  SimpleSample_load_intermediate(self$data_location)
}

# Corresponds to @staticmethod SimpleSample.load_intermediate: Read JSON, pickle or its ZIP.
SimpleSample_load_intermediate <- function(data_location) {
  if (length(data_location) != 1L || !is.character(data_location) ||
      !file.exists(data_location)) {
    stop("Failed to load sample data from: ", data_location, call. = FALSE)
  }

  sample_data <- NULL
  if (.samples_is_zipfile(data_location)) {
    first_entry <- .samples_read_first_zip_entry(data_location)
    entry_name <- first_entry$name

    # Python determines the decoder based on the extension of the first item in the ZIP.
    if (!is.null(entry_name) && endsWith(entry_name, ".pickle")) {
      sample_data <- .samples_read_pickle(data_location, entry_name)
    } else if (!is.null(entry_name) && endsWith(entry_name, ".json")) {
      json_text <- rawToChar(first_entry$raw)
      sample_data <- .samples_parse_json(json_text)
    }
  } else if (endsWith(data_location, ".pickle")) {
    sample_data <- .samples_read_pickle(data_location)
  } else if (endsWith(data_location, ".json")) {
    json_text <- paste(readLines(data_location, warn = FALSE), collapse = "\n")
    sample_data <- .samples_parse_json(json_text)
  } else {
    stop("Unknown file format: ", data_location, call. = FALSE)
  }

  if (is.null(sample_data)) {
    stop("Failed to load sample data from: ", data_location, call. = FALSE)
  }
  sample_data
}

# Use a constructor alias close to the Python class name for direct use by sample_factory in constructors.R.
SimpleSample <- SimpleSample__init__

# Preserves the compatibility constructor used by earlier smoke_check.R; it is not part of Python's 9 defs.
new_sample <- function(name, infile) {
  list(
    name = name,
    infile = infile,
    mass_tracks = NULL,
    peaks = NULL
  )
}
