# Corresponds to Python asari/db/__init__.py and the R-compatible positioning layer of its three data resources.

ASARI_DB_RESOURCES <- c(
  "gcms_libraries.json", "mass_indexed_compounds.pickle", "emp_cpds_trees.pickle"
)

# Find the data directory: first use the explicit option, then check the installed R package and the Python environment of the current project.
.asari_db_directories <- function() {
  explicit <- getOption("asariR.db_dir", Sys.getenv("ASARI_R_DB_DIR", ""))
  installed <- system.file("db", package = "asariR")
  python_root <- getOption(
    "asariR.python_package_root", Sys.getenv("ASARI_PYTHON_PACKAGE_ROOT", "")
  )
  python_explicit <- if (nzchar(python_root)) file.path(python_root, "db") else character()
  project_python <- c(
    Sys.glob(file.path(
      ".venv", "lib", "python*", "site-packages", "asari", "db"
    )),
    # Tests are typically executed from asari-r/, with the Python environment located in the parent project directory.
    Sys.glob(file.path(
      "..", ".venv", "lib", "python*", "site-packages", "asari", "db"
    ))
  )
  unique(Filter(nzchar, c(explicit, installed, python_explicit, project_python)))
}

# Returns the actual path to the specified database resource without using a hidden default data file.
asari_db_path <- function(name, mustWork = TRUE) {
  if (!name %in% ASARI_DB_RESOURCES) {
    stop("Unknown asari database resource: ", name, call. = FALSE)
  }
  directories <- .asari_db_directories()
  candidates <- file.path(directories, name)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) return(normalizePath(existing[[1L]]))
  if (isTRUE(mustWork)) {
    stop(
      "Database resource not found: ", name,
      ". Set options(asariR.db_dir = '/explicit/path/to/asari/db').",
      call. = FALSE
    )
  }
  if (length(candidates)) candidates[[1L]] else file.path("db", name)
}

# Read JSON or Python pickle in original format; pickle specifically converts and reuses the decoder of samples.R.
load_asari_db_resource <- function(name) {
  path <- asari_db_path(name)
  if (endsWith(name, ".json")) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Reading JSON database resources requires jsonlite.", call. = FALSE)
    }
    return(jsonlite::fromJSON(path, simplifyVector = FALSE))
  }
  if (!exists(".samples_read_pickle", mode = "function", inherits = TRUE)) {
    stop("Reading Python pickle database resources requires samples.R.", call. = FALSE)
  }
  .samples_read_pickle(path)
}
