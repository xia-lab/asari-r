# Corresponds to Python asari/tools/mzml_converter.py: Install and call ThermoRawFileParser in batches.

# Class-level version table and installation status.
mzMLconverter_versions <- list(
  `1.4.5` = "https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.5/ThermoRawFileParser1.4.5.zip"
)
.mzMLconverter_installed <- new.env(parent = emptyenv())

# Create a variable object corresponding to mzMLconverter.
mzMLconverter <- function(version = "1.4.5", dask_ip = FALSE, multicores = NULL) {
  mzMLconverter__init__(version, dask_ip, multicores)
}

# Corresponds to mzMLconverter.__init__.
mzMLconverter__init__ <- function(version = "1.4.5", dask_ip = FALSE, multicores = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("mzMLconverter", "environment")
  self$version <- version
  self$dask_ip <- dask_ip
  self$multicores <- if (is.null(multicores)) parallel::detectCores() else as.integer(multicores)
  mzMLconverter_install_converter(version)
  self$command_template <- mzMLconverter__conversion_command_template(self)
  self$bulk_convert <- function(raw_files) mzMLconverter_bulk_convert(self, raw_files)
  self
}

# Corresponds to internal nested def __determine_executable.
mzMLconverter__determine_executable <- function(converter_path) {
  probe <- getOption("asariR.converter_probe")
  engines <- list("mono", character())
  for (engine in engines) {
    ok <- if (is.function(probe)) {
      isTRUE(probe(c(engine, converter_path)))
    } else {
      executable <- if (length(engine)) Sys.which(engine[[1L]]) else converter_path
      nzchar(executable) && suppressWarnings(system2(executable, if (length(engine)) converter_path else character(), stdout = FALSE, stderr = FALSE)) %in% c(0L, 1L)
    }
    if (ok) return(engine)
  }
  stop("Could not determine the executable to run the converter!")
}

# Corresponds to the private method __conversion_command_template.
mzMLconverter__conversion_command_template <- function(self) {
  path <- if (exists(self$version, envir = .mzMLconverter_installed, inherits = FALSE)) {
    .mzMLconverter_installed[[self$version]]
  } else mzMLconverter_install_converter(self$version)
  converter <- file.path(path, "ThermoRawFileParser.exe")
  engine <- mzMLconverter__determine_executable(converter)
  c(engine, converter, "-i", "INPUT", "-o", "OUTPUT", "-f", "2")
}

# Corresponds to the static method install_converter.
mzMLconverter_install_converter <- function(version = "1.4.5") {
  if (exists(version, envir = .mzMLconverter_installed, inherits = FALSE)) {
    return(.mzMLconverter_installed[[version]])
  }
  url <- mzMLconverter_versions[[version]]
  if (is.null(url)) stop("Unsupported converter version: ", version)
  root <- getOption("asariR.converter_dir", tools::R_user_dir("asariR", "data"))
  extract_path <- file.path(root, paste0("ThermoRawFileParser", version))
  if (dir.exists(extract_path)) {
    .mzMLconverter_installed[[version]] <- normalizePath(extract_path)
    return(.mzMLconverter_installed[[version]])
  }
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(pattern = paste0("ThermoRawFileParser", version), fileext = ".zip", tmpdir = root)
  on.exit(if (file.exists(archive)) unlink(archive), add = TRUE)
  downloader <- getOption("asariR.converter_download")
  if (is.function(downloader)) downloader(url, archive)
  else utils::download.file(url, archive, mode = "wb", quiet = TRUE)
  utils::unzip(archive, exdir = extract_path)
  if (!dir.exists(extract_path)) stop("Failed to install the converter")
  .mzMLconverter_installed[[version]] <- normalizePath(extract_path)
  .mzMLconverter_installed[[version]]
}

# Corresponds to the static method uninstall_converter.
mzMLconverter_uninstall_converter <- function(version = "1.4.5") {
  if (exists(version, envir = .mzMLconverter_installed, inherits = FALSE)) {
    path <- .mzMLconverter_installed[[version]]
    if (dir.exists(path)) unlink(path, recursive = TRUE) else if (file.exists(path)) unlink(path)
    rm(list = version, envir = .mzMLconverter_installed)
  }
  invisible(NULL)
}

# Corresponds to bulk_convert: replaces INPUT/OUTPUT and executes system commands in parallel.
mzMLconverter_bulk_convert <- function(self, raw_files) {
  commands <- lapply(raw_files, function(file) {
    template <- gsub("INPUT", file, self$command_template, fixed = TRUE)
    gsub("OUTPUT", dirname(file), template, fixed = TRUE)
  })
  runner <- getOption("asariR.converter_run")
  if (!is.function(runner)) runner <- function(command) system2(command[[1L]], command[-1L])
  workers <- if (identical(self$dask_ip, FALSE)) self$multicores else self$dask_ip
  results <- bulk_process(runner, commands, jobs_per_worker = workers)
  for (ii in seq_along(results)) {
    if (identical(as.integer(results[[ii]]), 0L)) cat("Successfully converted ", raw_files[[ii]], "\n", sep = "")
    else cat("Failed to convert ", raw_files[[ii]], "\n", sep = "")
  }
  invisible(NULL)
}
