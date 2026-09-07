# Corresponds to Python asari/tools/gc_annotation.py: old version of GC annotation class and EI spectral library package.

# Create a GC_Annotation mix-in object; Python classes do not have __init__, and the R factory is only responsible for providing mutable containers.
GC_Annotation <- function(parameters = list(), sample_registry = list(), sample_factory = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("GC_Annotation", "environment")
  self$parameters <- parameters
  self$sample_registry <- sample_registry
  self$sample_factory <- sample_factory
  self$RI_map <- NULL
  self$RI_models <- list()
  self$reverse_RI_models <- list()
  self$populate_RI_lookup <- function(sample_map) GC_Annotation_populate_RI_lookup(self, sample_map)
  self$convert_to_RI <- function(sample_map) GC_Annotation_convert_to_RI(self, sample_map)
  self$annotate_GC <- function() GC_Annotation_annotate_GC(self)
  self
}

# Construct an interpolation function that supports linear extrapolation outside the interval, corresponding to scipy interp1d(..., extrapolate).
.gc_approxfun_extrapolate <- function(x, y) {
  ordering <- order(x)
  x <- as.numeric(x[ordering])
  y <- as.numeric(y[ordering])
  function(newx) {
    result <- stats::approx(x, y, xout = newx, rule = 1, ties = mean)$y
    if (length(x) >= 2L) {
      low <- newx < x[[1L]]
      high <- newx > x[[length(x)]]
      result[low] <- y[[1L]] + (newx[low] - x[[1L]]) * (y[[2L]] - y[[1L]]) / (x[[2L]] - x[[1L]])
      last <- length(x)
      result[high] <- y[[last]] + (newx[high] - x[[last]]) *
        (y[[last]] - y[[last - 1L]]) / (x[[last]] - x[[last - 1L]])
    }
    result
  }
}

# Read fields from environment or list samples.
.gc_sample_field <- function(sample, names) {
  for (name in names) {
    value <- if (is.environment(sample)) sample[[name]] else sample[[name]]
    if (!is.null(value)) return(value)
  }
  NULL
}

# Corresponds to GC_Annotation.populate_RI_lookup: Fit RT to RI and RI to scan number model from standard RT.
GC_Annotation_populate_RI_lookup <- function(self, sample_map) {
  ri_maps <- list()
  ri_models <- list()
  reverse_models <- list()
  standards <- utils::read.csv(self$parameters$retention_index_standards, check.names = FALSE)
  reference_ids <- unique(unname(unlist(sample_map, use.names = FALSE)))
  factory <- self$sample_factory
  if (!is.function(factory)) factory <- getOption("asariR.gc_sample_factory")
  if (!is.function(factory)) stop("GC RI calibration requires a sample factory.")

  for (reference_id in reference_ids) {
    key <- as.character(reference_id)
    cat(key, "\n")
    ri_maps[[key]] <- list()
    sample <- factory(self$sample_registry[[key]], experiment = self)
    sample_name <- .gc_sample_field(sample, "name")
    retention_times <- as.numeric(.gc_sample_field(sample, c("list_retention_time", "retention_time")))
    scan_numbers <- as.numeric(.gc_sample_field(sample, c("list_scan_numbers", "rt_numbers")))
    previous_index <- NULL
    next_index <- NULL
    previous_rt <- NULL
    next_rt <- NULL
    indexes <- numeric()
    for (ii in seq_along(retention_times)) {
      rt <- retention_times[[ii]]
      scan_no <- scan_numbers[[ii]]
      cat(rt, scan_no, "\n")
      for (jj in seq_len(nrow(standards))) {
        index <- as.integer(standards$Index[[jj]])
        index_rt <- as.numeric(standards[[sample_name]][[jj]])
        cat("\t", index, index_rt, "\n")
        if (rt > index_rt) {
          previous_index <- index
          previous_rt <- index_rt
        } else {
          next_index <- index
          next_rt <- index_rt
          break
        }
      }
      if (is.null(next_rt)) {
        next_rt <- max(retention_times) * 1.1
        next_index <- max(standards$Index) + 1
      }
      invisible(next_index)
      if (is.null(previous_index) || is.null(previous_rt)) {
        stop("Retention time precedes the first RI standard, matching Python failure.")
      }
      # Keep the (next_rt - rt) denominator in the original Python formula.
      ri_value <- 100 * (previous_index + ((rt - previous_rt) / (next_rt - rt)))
      indexes <- c(indexes, ri_value)
      ri_maps[[key]][[as.character(rt)]] <- ri_value
    }
    model <- stats::lowess(retention_times, indexes)
    model2 <- stats::lowess(scan_numbers, indexes)
    ri_models[[key]] <- .gc_approxfun_extrapolate(model$x, model$y)
    reverse_models[[key]] <- .gc_approxfun_extrapolate(model2$y, model2$x)
  }
  self$RI_map <- ri_maps
  self$RI_models <- ri_models
  self$reverse_RI_models <- reverse_models
  invisible(NULL)
}

# Corresponds to GC_Annotation.convert_to_RI: write to the registry according to the model referenced by sample_map.
GC_Annotation_convert_to_RI <- function(self, sample_map) {
  if (is.null(self$RI_map) || length(self$RI_map) == 0L) self$populate_RI_lookup(sample_map)
  for (key in names(sample_map)) {
    sample <- self$sample_registry[[key]]
    retention_times <- sample$list_retention_time
    sample$list_retention_index <- self$RI_models[[as.character(sample_map[[key]])]](retention_times)
    self$sample_registry[[key]] <- sample
  }
  invisible(NULL)
}

# Corresponds to GC_Annotation.annotate_GC: annotate the preferred and full tables.
GC_Annotation_annotate_GC <- function(self) {
  preferred <- file.path(
    self$parameters$outdir,
    paste0("preferred_", self$parameters$output_feature_table)
  )
  full <- file.path(
    self$parameters$outdir, "export",
    paste0("full_", self$parameters$output_feature_table)
  )
  EI_MS_Library_annotate_gc_feature_table_with_library(preferred, self$parameters$GC_Database)
  EI_MS_Library_annotate_gc_feature_table_with_library(full, self$parameters$GC_Database)
  invisible(NULL)
}

# Create a variable object corresponding to EI_MS_Library.
EI_MS_Library <- function(library_ID, multicores = NULL) EI_MS_Library__init__(library_ID, multicores)

# Corresponds to EI_MS_Library.__init__.
EI_MS_Library__init__ <- function(library_ID, multicores = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("EI_MS_Library", "environment")
  self$library_meta <- EI_MS_Library_retrieve_library_meta(self, library_ID)
  self$library <- NULL
  cores <- parallel::detectCores()
  self$multicores <- if (is.null(multicores)) cores else min(as.integer(multicores), cores)
  self$load_library <- function(limit = NULL) EI_MS_Library_load_library(self, limit)
  self$retrieve_library_meta <- function(id) EI_MS_Library_retrieve_library_meta(self, id)
  self$annotate_gc_feature_table <- function(feature_table_path, drt = 0.5, min_peaks = 3L,
                                              min_shared_peaks = 1L, min_score_threshold = 0.7) {
    EI_MS_Library_annotate_gc_feature_table(
      self, feature_table_path, drt, min_peaks, min_shared_peaks, min_score_threshold
    )
  }
  self$load_library()
  self
}

# Unify MSP/MGF/JSON records into lightweight Spectrum objects.
.gc_library_spectrum <- function(entry) {
  if (!is.null(entry$mz) && !is.null(entry$intensities)) return(entry)
  peaks <- entry$peaks
  matrix_peaks <- if (is.matrix(peaks)) peaks else if (length(peaks)) {
    do.call(rbind, lapply(peaks, function(peak) {
      if (is.list(peak) && !is.null(peak$mz)) c(peak$mz, peak$intensity) else as.numeric(peak[1:2])
    }))
  } else matrix(numeric(), ncol = 2L)
  metadata <- entry
  metadata$peaks <- NULL
  .feature_graph_spectrum(matrix_peaks[, 1L], matrix_peaks[, 2L], metadata)
}

# Corresponds to EI_MS_Library.load_library: select parser, limit number and filter normalization.
EI_MS_Library_load_library <- function(self, limit = NULL) {
  extension <- sub("^\\.", "", self$library_meta$Extension)
  loader <- getOption("asariR.spectral_library_loader")
  if (is.function(loader)) {
    library <- loader(self$library_meta$LIB_PATH, extension)
  } else if (tolower(extension) == "msp") {
    library <- parse_msp(self$library_meta$LIB_PATH)
  } else if (tolower(extension) == "mgf") {
    library <- parse_mgf(self$library_meta$LIB_PATH)
  } else if (tolower(extension) == "json" && requireNamespace("jsonlite", quietly = TRUE)) {
    library <- jsonlite::fromJSON(self$library_meta$LIB_PATH, simplifyVector = FALSE)
  } else {
    stop("No MatchMS-compatible loader for extension: ", extension)
  }
  if (!is.null(limit)) library <- utils::head(library, as.integer(limit) + 1L)
  self$loader <- loader
  self$library <- lapply(library, function(entry) {
    .feature_graph_process_spectrum(.gc_library_spectrum(entry))
  })
  invisible(NULL)
}

# Locate the asari database directory.
.gc_annotation_db_dir <- function() {
  configured <- getOption("asariR.db_dir", "")
  if (nzchar(configured)) return(configured)
  installed <- system.file("db", package = "asariR")
  if (nzchar(installed)) return(installed)
  candidates <- c(
    file.path(".venv", "lib", "python3.13", "site-packages", "asari", "db"),
    file.path("..", ".venv", "lib", "python3.13", "site-packages", "asari", "db")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing)) existing[[1L]] else ""
}

# Corresponds to EI_MS_Library.retrieve_library_meta.
EI_MS_Library_retrieve_library_meta <- function(self, library_ID) {
  invisible(self)
  if (file.exists(library_ID) && !dir.exists(library_ID)) {
    cat("Assuming provided file is the library\n")
    cat("Assuming MatchMS is okay for parsing if in JSON, MSP, or MGF format\n")
    return(list(
      Name = library_ID,
      Description = "User Provided Library",
      URL = NULL,
      OnDiskName = basename(library_ID),
      Extension = paste0(".", tools::file_ext(library_ID)),
      Parser = "matchms",
      License = "Unknown",
      `Required Citations` = "Unknown",
      LIB_PATH = library_ID
    ))
  }
  manifest <- EI_MS_Library_load_library_manifest()
  entry <- manifest$EI_MS[[library_ID]]
  if (is.null(entry)) {
    cat("Library ID ", library_ID, " not found\n", sep = "")
    cat("Valid Selections:\n\t", paste(names(manifest$EI_MS), collapse = "\n\t"), "\n", sep = "")
    stop("Unknown library ID: ", library_ID)
  }
  on_disk <- file.path(.gc_annotation_db_dir(), entry$OnDiskName)
  if (!file.exists(on_disk)) {
    cat("Downloading Library...\n")
    if (!is.null(entry$URL) && (endsWith(entry$URL, "zip") || identical(entry$Compression, "zip"))) {
      download_and_unzip_to_pkg_resources(entry$URL, "asariR", "db")
    }
    if (!file.exists(on_disk)) stop("Library not found at ", on_disk, " after Download")
  }
  if (!is.null(entry$LIB_PATH)) stop("Library already loaded!")
  entry$LIB_PATH <- on_disk
  entry
}

# Corresponds to the static method load_library_manifest; option is used for test injection, and the disk results are cached in the environment.
.gc_manifest_cache <- new.env(parent = emptyenv())
EI_MS_Library_load_library_manifest <- function() {
  injected <- getOption("asariR.gcms_library_manifest")
  if (!is.null(injected)) return(injected)
  if (exists("manifest", envir = .gc_manifest_cache, inherits = FALSE)) {
    return(.gc_manifest_cache$manifest)
  }
  path <- file.path(.gc_annotation_db_dir(), "gcms_libraries.json")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading GC library manifest requires jsonlite.")
  .gc_manifest_cache$manifest <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  .gc_manifest_cache$manifest
}

# Corresponds to annotate_gc_feature_table: Construct co-elution spectra and compare them one by one with the entire EI library.
EI_MS_Library_annotate_gc_feature_table <- function(
    self,
    feature_table_path,
    drt = 0.5,
    min_peaks = 3L,
    min_shared_peaks = 1L,
    min_score_threshold = 0.7) {
  raw <- FeatureGraph_ftgraph_from_ft(feature_table_path)
  coelute <- raw$filter_graph(drt)
  extracted <- coelute$extract_fragmentation_spectrum(TRUE, min_peaks)
  cat("Total Cluster Spectra: ", length(extracted), "\n", sep = "")
  cat("Total Library Spectra: ", length(self$library), "\n", sep = "")
  cat("Total Comparisons: ", length(extracted) * length(self$library), ", this may take some time...\n", sep = "")
  matches <- list()
  for (extract in extracted) {
    for (library in self$library) {
      result <- wrapped_cosine(list(extract, library))
      score <- as.numeric(result[[2L]])
      if (length(score) == 2L && score[[1L]] >= min_score_threshold && score[[2L]] >= min_shared_peaks) {
        matches[[length(matches) + 1L]] <- list(
          extract = extract, library = library,
          similarity = score[[1L]], match_peaks = score[[2L]]
        )
      }
    }
  }
  coelute$map_annotations(matches)
  invisible(NULL)
}

# Corresponds to the static method annotate_gc_feature_table_with_library.
EI_MS_Library_annotate_gc_feature_table_with_library <- function(
    feature_table_path, library_ID, multicores = NULL) {
  library <- EI_MS_Library(library_ID, multicores)
  library$annotate_gc_feature_table(feature_table_path)
  invisible(NULL)
}

# Corresponds to wrapped_cosine: Returns the original job and CosineGreedy score tuples.
wrapped_cosine <- function(job) {
  pairer <- getOption("asariR.cosine_greedy_pair")
  score <- if (is.function(pairer)) {
    pairer(job[[1L]], job[[2L]])
  } else {
    cosine_similarity(job[[1L]]$peaks, job[[2L]]$peaks)
  }
  list(job, unlist(score, use.names = FALSE))
}
