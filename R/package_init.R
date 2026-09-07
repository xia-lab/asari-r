# Corresponds to Python asari/__init__.py: declare the package version and load each module uniformly by the R package loading mechanism.

# Python provides externally __version__ = "1.17.0"; backticks allow R to retain the symbol with the same name.
`__version__` <- "1.17.0"

# Provide aliases with the same value that conform to R naming conventions.
asari_version <- `__version__`

# The import of Python __init__ will load the main module; the R package will automatically read all files under R/ when loading.
ASARI_CORE_MODULES <- c(
  "analyze", "annotate", "chromatograms", "command_line", "constructors",
  "default_parameters", "experiment", "main", "mass_functions", "peaks",
  "samples", "workflow"
)
