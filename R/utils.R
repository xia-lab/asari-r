# Corresponds to Python asari/utils.py: parallel tasks, file downloads, checksums, and mzML tools.

# Convert an R matrix to a nested list in row order with Python ndarray.tolist().
.utils_matrix_to_list <- function(value) {
  lapply(seq_len(nrow(value)), function(row) {
    as.list(unname(value[row, , drop = TRUE]))
  })
}

# Download URL to temporary file; file:// branch is only for local verification without network.
.utils_download <- function(url, destination) {
  if (startsWith(url, "file://")) {
    source <- sub("^file://", "", url)
    if (!file.copy(source, destination, overwrite = TRUE)) {
      stop("Failed to copy local archive: ", source, call. = FALSE)
    }
    return(invisible(destination))
  }
  status <- utils::download.file(
    url, destination, mode = "wb", quiet = TRUE
  )
  if (!identical(status, 0L)) {
    stop("Failed to download archive: ", url, call. = FALSE)
  }
  invisible(destination)
}

# Locate the R package resource directory; directory path input is used for local development and testing.
.utils_package_path <- function(package) {
  if (dir.exists(package)) return(normalizePath(package))
  path <- system.file(package = package)
  if (!nzchar(path)) stop("R package is not installed: ", package, call. = FALSE)
  path
}

# Get mzML scan header information from the path or header provided by the test.
.utils_mzml_header <- function(mzml_file) {
  if (is.list(mzml_file) && !is.null(mzml_file$header)) {
    return(mzml_file$header)
  }
  if (is.environment(mzml_file) && exists(
      "header", envir = mzml_file, inherits = FALSE)) {
    return(mzml_file$header)
  }
  path <- if (is.character(mzml_file)) mzml_file else {
    if (is.environment(mzml_file)) mzml_file$path else mzml_file$path
  }
  handle <- mzR::openMSfile(path)
  on.exit(mzR::close(handle), add = TRUE)
  mzR::header(handle)
}

# Windows does not have fork, and the temporary PSOCK cluster is used alone to perform parallel tasks.
.utils_psock_apply <- function(arguments, command, workers) {
  worker_cluster <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(worker_cluster), add = TRUE)
  parallel::parLapply(worker_cluster, arguments, command)
}

# Creates a lightweight mutable object corresponding to a Python json.JSONEncoder subclass.
NpEncoder <- function() {
  self <- new.env(parent = emptyenv())
  class(self) <- c("NpEncoder", "environment")
  self$default <- function(obj) NpEncoder_default(self, obj)
  self
}

# Corresponds to NpEncoder.default: Convert the numpy corresponding type into a JSON-encodable basic type.
NpEncoder_default <- function(self, obj) {
  # self corresponds to a Python JSONEncoder instance, and the conversion itself does not rely on instance state.
  invisible(self)
  if (is.integer(obj) && length(obj) == 1L) return(as.integer(obj))
  if (is.numeric(obj) && length(obj) == 1L) return(as.numeric(obj))
  if (is.matrix(obj)) return(.utils_matrix_to_list(obj))
  if (is.array(obj) || (is.atomic(obj) && length(obj) > 1L)) {
    return(as.list(unname(obj)))
  }
  stop(
    "Object of class ", paste(class(obj), collapse = "/"),
    " is not JSON serializable.",
    call. = FALSE
  )
}

# Corresponds to bulk_process: applies command to arguments one by one and returns the results in the order of input.
bulk_process <- function(command, arguments, jobs_per_worker = FALSE) {
  if (!is.function(command)) stop("command must be a function.", call. = FALSE)
  if (!is.list(arguments)) arguments <- as.list(arguments)
  workers <- if (
      isFALSE(jobs_per_worker) || identical(jobs_per_worker, "auto")) {
    parallel::detectCores()
  } else {
    as.integer(jobs_per_worker)
  }
  if (is.na(workers) || workers < 1L) stop("Invalid worker count.")
  if (length(arguments) == 0L) return(list())
  if (workers == 1L) return(lapply(arguments, command))

  # Unix uses fork; Windows uses temporary PSOCK cluster, both retaining the input order of starmap.
  if (.Platform$OS.type != "windows") {
    parallel::mclapply(arguments, command, mc.cores = workers)
  } else {
    .utils_psock_apply(arguments, command, workers)
  }
}

# Corresponds to download_and_unzip_to_pkg_resources: Download ZIP to a subdirectory of the package parent directory.
download_and_unzip_to_pkg_resources <- function(url, package, subdir = "data") {
  message("HERE")
  package_dir <- dirname(.utils_package_path(package))
  extract_to <- file.path(package_dir, subdir)
  dir.create(extract_to, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".zip")
  on.exit(unlink(archive), add = TRUE)
  .utils_download(url, archive)
  utils::unzip(archive, exdir = extract_to)
  message("Extracted to: ", extract_to)
  invisible(NULL)
}

# Corresponds to download_and_unzip: download the ZIP and extract it to the specified target directory.
download_and_unzip <- function(url, extract_to) {
  dir.create(extract_to, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".zip")
  on.exit(unlink(archive), add = TRUE)
  .utils_download(url, archive)
  utils::unzip(archive, exdir = extract_to)
  message("Extracted to: ", extract_to)
  invisible(NULL)
}

# Corresponds to validate_mzml_file: complete mzML scan, any read exception returns FALSE.
validate_mzml_file <- function(file) {
  tryCatch({
    handle <- mzR::openMSfile(file)
    on.exit(mzR::close(handle), add = TRUE)
    scan_header <- mzR::header(handle)
    if (nrow(scan_header) > 0L) {
      for (scan in seq_len(nrow(scan_header))) mzR::peaks(handle, scan)
    }
    TRUE
  }, error = function(error) FALSE)
}

# Corresponds to build_boolean_dict: holds the 8 Boolean input keys accepted by Python.
build_boolean_dict <- function() {
  list(
    T = TRUE,
    F = FALSE,
    `1` = TRUE,
    `0` = FALSE,
    True = TRUE,
    False = FALSE,
    `TRUE` = TRUE,
    `FALSE` = FALSE,
    true = TRUE,
    false = FALSE
  )
}

# Corresponds to sizeof_fmt: Format a number of bytes into text in the IEC units used by Python.
sizeof_fmt <- function(num, suffix = "B") {
  units <- c("", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi")
  value <- as.numeric(num)
  for (unit in units) {
    if (abs(value) < 1024) return(sprintf("%3.1f%s%s", value, unit, suffix))
    value <- value / 1024
  }
  sprintf("%.1fYi%s", value, suffix)
}

# Corresponds to checksum_file: calculate file MD5; chunksize retains Python interface.
checksum_file <- function(file, chunksize = 16384) {
  if (!file.exists(file) || dir.exists(file)) stop("file must be a regular file.")
  if (length(chunksize) != 1L || chunksize <= 0) stop("Invalid chunksize.")
  unname(tools::md5sum(file))
}

# Corresponds to wait_with_pbar: the waiting progress bar is updated every second.
wait_with_pbar <- function(wait = 5) {
  if (length(wait) != 1L || is.na(wait) || wait < 0 || wait != as.integer(wait)) {
    stop("wait must be a non-negative integer.")
  }
  wait <- as.integer(wait)
  if (wait == 0L) return(invisible(NULL))
  progress <- utils::txtProgressBar(min = 0, max = wait, style = 3)
  on.exit(close(progress), add = TRUE)
  for (ii in seq_len(wait)) {
    Sys.sleep(1)
    utils::setTxtProgressBar(progress, ii)
  }
  invisible(NULL)
}

# Corresponds to get_ionization_mode_mzml: Check the positive and negative ion mode of limit+1 scans before.
get_ionization_mode_mzml <- function(mzml_file, limit = 50) {
  header <- .utils_mzml_header(mzml_file)
  if (!("polarity" %in% names(header))) stop("mzML polarity is unavailable.")
  modes <- character()
  scan_limit <- min(nrow(header), as.integer(limit) + 1L)
  if (scan_limit > 0L) {
    for (ii in seq_len(scan_limit)) {
      polarity <- header$polarity[[ii]]
      if (is.na(polarity) || !(polarity %in% c(0, 1))) {
        stop("Scan has neither positive nor negative polarity.")
      }
      modes <- unique(c(modes, if (polarity == 1) "pos" else "neg"))
      if (length(modes) > 1L) return("mixed")
    }
  }
  if (length(modes) == 0L) stop("mzML contains no spectra.")
  modes[[1L]]
}

# The ppm conversion helper function used by earlier R modules is retained; it is not part of Python's 10 defs.
ppm_to_mz_tolerance <- function(mz, ppm) {
  mz * ppm * 1e-6
}
