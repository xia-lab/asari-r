# Corresponds to Python asari/tools/file_io.py: Conversion between feature table and JSON style records.

# Corresponds to read_features_from_asari_table: reads the peak list from the asari table text.
read_features_from_asari_table <- function(
    text_table,
    id_col = 0L,
    mz_col = 1L,
    rtime_col = 2L,
    left_base = 3L,
    right_base = 4L,
    parent_masstrack_id = 5L,
    peak_area = 6L,
    cSelectivity = 7L,
    goodness_fitting = 8L,
    snr = 9L,
    detection_counts = 10L,
    delimiter = "\t") {
  # The Python parameter is the base number of 0, which is increased by 1 after entering R.
  column <- as.integer(c(
    id_col, mz_col, rtime_col, left_base, right_base,
    parent_masstrack_id, peak_area, cSelectivity,
    goodness_fitting, snr, detection_counts
  )) + 1L
  feature_lines <- strsplit(text_table, "\n", fixed = TRUE)[[1L]]
  if (length(feature_lines) == 0L || !nzchar(feature_lines[[1L]])) {
    stop("Feature table is empty.", call. = FALSE)
  }
  split_fixed <- function(line) strsplit(line, delimiter, fixed = TRUE)[[1L]]
  header <- split_fixed(sub("\r$", "", feature_lines[[1L]]))
  num_features <- length(feature_lines) - 1L
  num_samples <- length(header) - 11L
  cat("table header looks like: \n  ", paste(utils::head(header, 20L), collapse = ", "), "\n")
  cat(sprintf("Read %d feature lines\n", num_features))

  peaks <- list()
  if (num_features > 0L) {
    for (ii in seq.int(2L, length(feature_lines))) {
      line <- sub("\r$", "", feature_lines[[ii]])
      if (!nzchar(trimws(line))) next
      fields <- split_fixed(line)
      if (length(fields) < max(column)) stop("Feature line has too few columns.")
      feature_id <- fields[[column[[1L]]]]
      peaks[[length(peaks) + 1L]] <- list(
        id_number = feature_id,
        id = feature_id,
        mz = as.numeric(fields[[column[[2L]]]]),
        rtime = as.numeric(fields[[column[[3L]]]]),
        apex = as.numeric(fields[[column[[3L]]]]),
        rtime_left_base = as.numeric(fields[[column[[4L]]]]),
        rtime_right_base = as.numeric(fields[[column[[5L]]]]),
        parent_masstrack_id = fields[[column[[6L]]]],
        peak_area = as.numeric(fields[[column[[7L]]]]),
        cSelectivity = as.numeric(fields[[column[[8L]]]]),
        goodness_fitting = as.numeric(fields[[column[[9L]]]]),
        snr = as.numeric(fields[[column[[10L]]]]),
        detection_counts = as.integer(fields[[column[[11L]]]])
      )
    }
  }
  list(num_samples, peaks)
}

# Corresponds to export_json_to_table: the field name is used as the header, and the id field is always placed in the first column.
export_json_to_table <- function(j, outfile, sep = "\t") {
  if (length(j) == 0L) stop("j must contain at least one record.", call. = FALSE)
  fields <- names(j[[1L]])
  if ("id" %in% fields) fields <- c("id", fields[fields != "id"])
  lines <- c(
    paste(fields, collapse = sep),
    vapply(j, function(record) {
      paste(vapply(fields, function(field) as.character(record[[field]]), ""), collapse = sep)
    }, "")
  )
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(NULL)
}

# Corresponds to read_table_to_json: all fields are read as strings, and types are not automatically inferred.
read_table_to_json <- function(file, sep = "\t") {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) == 0L) stop("Table is empty.", call. = FALSE)
  header <- strsplit(sub("\r$", "", lines[[1L]]), sep, fixed = TRUE)[[1L]]
  if (length(lines) == 1L) return(list())
  lapply(lines[-1L], function(line) {
    values <- strsplit(sub("\r$", "", line), sep, fixed = TRUE)[[1L]]
    # Python zip will truncate the longer side; in R, only the common length of both is retained.
    size <- min(length(header), length(values))
    stats::setNames(as.list(values[seq_len(size)]), header[seq_len(size)])
  })
}
