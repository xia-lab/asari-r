# Corresponds to Python asari/command_line.py: This file only introduces the same entry from the main module.

# Use runtime lookup to avoid main.R not being read when the R package is loaded by file name.
command_line_main <- function(...) {
  entry <- get("main", mode = "function", inherits = TRUE)
  entry(...)
}
