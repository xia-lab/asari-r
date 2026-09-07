# Public entry points for table merging, GC feature graphs, and plotting.

# Identify the sample intensity column from the standard asari feature table.
.asari_sample_columns <- function(table) {
  metadata <- c(
    "id_number", "mz", "rtime", "rtime_left_base", "rtime_right_base",
    "parent_masstrack_id", "peak_area", "cSelectivity", "goodness_fitting",
    "snr", "detection_counts"
  )
  setdiff(names(table), metadata)
}

# Find the nearest available group within the current consensus group for a new feature.
.asari_merge_candidate <- function(clusters, mz, rtime, table_index, ppm, rt_tolerance) {
  if (length(clusters) == 0L) return(NA_integer_)
  centers_mz <- vapply(clusters, function(cluster) stats::median(cluster$mz), 0)
  delta_da <- mz * ppm * 1e-6
  ordering <- order(centers_mz)
  sorted <- centers_mz[ordering]
  left <- findInterval(mz - delta_da, sorted) + 1L
  right <- findInterval(mz + delta_da, sorted)
  if (left > right || left > length(sorted) || right < 1L) return(NA_integer_)
  candidates <- ordering[seq.int(max(1L, left), min(length(sorted), right))]
  candidates <- candidates[vapply(candidates, function(index) {
    cluster <- clusters[[index]]
    !(table_index %in% cluster$table_index) &&
      abs(rtime - stats::median(cluster$rtime)) <= rt_tolerance
  }, FALSE)]
  if (length(candidates) == 0L) return(NA_integer_)
  distance <- vapply(candidates, function(index) {
    cluster <- clusters[[index]]
    mass_part <- abs(mz - stats::median(cluster$mz)) / max(delta_da, .Machine$double.eps)
    rt_part <- abs(rtime - stats::median(cluster$rtime)) / rt_tolerance
    mass_part + rt_part
  }, 0)
  candidates[[which.min(distance)]]
}

#' Merge multiple asari feature tables
#'
#' `method = "standard"` directly merges ordinary asari tables; `method = "legacy"`
#' Call the old formula-mass table combiner corresponding to Python `tools/merge.py`.
#'
#' @param feature_tables Two or more explicit TSV paths.
#' @param output Explicit output path for the merged TSV file.
#' @param ppm m/z tolerance.
#' @param rt_tolerance RT tolerance, unit is seconds.
#' @param method Standard asari tables or old formula-mass tables.
#' @param make_feature_id Passed to legacy combiners only.
#' @return `standard` returns the merged data.frame; `legacy` returns the output path.
#' @export
asari_merge_feature_tables <- function(
    feature_tables,
    output,
    ppm = 4,
    rt_tolerance = 15,
    method = c("standard", "legacy"),
    make_feature_id = FALSE) {
  method <- match.arg(method)
  ppm <- .asari_public_positive_number(ppm, "ppm")
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  if (!is.character(feature_tables) || length(feature_tables) < 2L) {
    stop("feature_tables must contain at least two explicit TSV paths.", call. = FALSE)
  }
  paths <- vapply(feature_tables, function(path) {
    path <- path.expand(path)
    if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
    normalizePath(path, mustWork = TRUE)
  }, "")
  outfile <- .asari_public_output_file(output)
  if (method == "legacy") {
    merge_feature_tables(paths, make_feature_id, ppm, rt_tolerance, outfile)
    return(normalizePath(outfile, mustWork = TRUE))
  }

  tables <- lapply(paths, .asari_public_feature_table)
  invisible(lapply(seq_along(tables), function(ii) {
    .asari_public_require_columns(tables[[ii]], c("mz", "rtime"), paste0("feature_tables[[", ii, "]]"))
  }))
  for (ii in seq_along(tables)) {
    if (!("id_number" %in% names(tables[[ii]]))) {
      tables[[ii]]$id_number <- paste0("row", seq_len(nrow(tables[[ii]])))
    }
  }

  # Features are added to the consensus group in the order of the input tables, with each table occupying at most one row in a group.
  clusters <- list()
  for (table_index in seq_along(tables)) {
    table <- tables[[table_index]]
    row_order <- order(as.numeric(table$mz), as.numeric(table$rtime))
    for (row_index in row_order) {
      mz <- as.numeric(table$mz[[row_index]])
      rtime <- as.numeric(table$rtime[[row_index]])
      if (!is.finite(mz) || !is.finite(rtime)) next
      candidate <- .asari_merge_candidate(
        clusters, mz, rtime, table_index, ppm, rt_tolerance
      )
      member <- list(table_index = table_index, row_index = row_index)
      if (is.na(candidate)) {
        clusters[[length(clusters) + 1L]] <- list(
          mz = mz, rtime = rtime, table_index = table_index, members = list(member)
        )
      } else {
        clusters[[candidate]]$mz <- c(clusters[[candidate]]$mz, mz)
        clusters[[candidate]]$rtime <- c(clusters[[candidate]]$rtime, rtime)
        clusters[[candidate]]$table_index <- c(clusters[[candidate]]$table_index, table_index)
        clusters[[candidate]]$members <- c(clusters[[candidate]]$members, list(member))
      }
    }
  }
  clusters <- clusters[order(
    vapply(clusters, function(cluster) stats::median(cluster$mz), 0),
    vapply(clusters, function(cluster) stats::median(cluster$rtime), 0)
  )]

  # All sample columns are prefixed with the source table to avoid overwriting when different projects use the same sample name.
  table_labels <- make.unique(tools::file_path_sans_ext(basename(paths)))
  sample_names <- unlist(lapply(seq_along(tables), function(ii) {
    paste0(table_labels[[ii]], "::", .asari_sample_columns(tables[[ii]]))
  }), use.names = FALSE)
  rows <- lapply(seq_along(clusters), function(ii) {
    cluster <- clusters[[ii]]
    intensities <- stats::setNames(as.list(rep(NA_real_, length(sample_names))), sample_names)
    source_ids <- character()
    for (member in cluster$members) {
      table_index <- member$table_index
      row_index <- member$row_index
      table <- tables[[table_index]]
      source_ids <- c(source_ids, paste0(
        table_labels[[table_index]], ":", as.character(table$id_number[[row_index]])
      ))
      for (sample in .asari_sample_columns(table)) {
        key <- paste0(table_labels[[table_index]], "::", sample)
        intensities[[key]] <- suppressWarnings(as.numeric(table[[sample]][[row_index]]))
      }
    }
    c(list(
      consensus_id = paste0("C", ii - 1L),
      mz = stats::median(cluster$mz), rtime = stats::median(cluster$rtime),
      detection_tables = length(unique(cluster$table_index)),
      source_ids = paste(source_ids, collapse = ";")
    ), intensities)
  })
  result <- if (length(rows)) {
    as.data.frame(do.call(rbind, lapply(rows, function(row) {
      as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
    })), check.names = FALSE)
  } else data.frame(
    consensus_id = character(), mz = numeric(), rtime = numeric(),
    detection_tables = integer(), source_ids = character()
  )
  # rbind will convert the mixed list into text, and restore the consensus numerical column here.
  for (name in c("mz", "rtime", "detection_tables", sample_names)) {
    if (name %in% names(result)) result[[name]] <- suppressWarnings(as.numeric(result[[name]]))
  }
  utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  attr(result, "output_file") <- outfile
  result
}

#' Construct GC feature graph according to co-elution relationship
#'
#' @param feature_table an explicit asari feature table TSV.
#' @param rt_tolerance Build the RT tolerance for edges.
#' @param output Optional feature - Cluster TSV output path.
#' @return A data.frame with one feature per row and its cluster ID.
#' @export
asari_feature_graph <- function(feature_table, rt_tolerance = 0.5, output = NULL) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  rt_tolerance <- .asari_public_positive_number(rt_tolerance, "rt_tolerance")
  graph <- FeatureGraph_ftgraph_from_ft(normalizePath(path, mustWork = TRUE))$filter_graph(rt_tolerance)
  graph$find_spectral_clusters()
  rows <- unlist(lapply(names(graph$clusters), function(cluster_id) {
    lapply(graph$clusters[[cluster_id]], function(feature_id) {
      data.frame(feature_id = feature_id, cluster_id = as.integer(cluster_id))
    })
  }), recursive = FALSE)
  result <- if (length(rows)) do.call(rbind, rows) else {
    data.frame(feature_id = character(), cluster_id = integer())
  }
  if (!is.null(output)) {
    outfile <- .asari_public_output_file(output)
    utils::write.table(result, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
    attr(result, "output_file") <- outfile
  }
  result
}

#' Draw a mirror image of two spectra
#'
#' @param query,reference Two columns of `m/z, intensity` matrix or data.frame.
#' @param output explicit PDF output path.
#' @param tolerance Optional fragment m/z matching tolerance Da.
#' @param title Figure title.
#' @return PDF absolute path.
#' @export
asari_mirror_plot <- function(query, reference, output, tolerance = NULL, title = "MS Mirror Plot") {
  outfile <- .asari_public_output_file(output)
  mirror_plot(
    as.matrix(query), as.matrix(reference), match_tol = tolerance,
    title = title, outfile = outfile
  )
  normalizePath(outfile, mustWork = TRUE)
}

#' Draw mzML specifying scan number and m/z area
#'
#' @param input An explicit mzML file.
#' @param output explicit PDF output path.
#' @param scan_range scan number open interval with length 2.
#' @param mz_range m/z open interval of length 2.
#' @param ms_level mass spectrum level.
#' @param style `"scatter"`, `"double"` or `"line"`.
#' @return PDF absolute path.
#' @export
asari_plot_mz_region <- function(
    input,
    output,
    scan_range,
    mz_range,
    ms_level = 1L,
    style = c("scatter", "double", "line")) {
  infile <- .asari_public_one_mzml(input)
  if (!is.numeric(scan_range) || length(scan_range) != 2L || diff(scan_range) <= 0) {
    stop("scan_range must be two increasing numbers.")
  }
  if (!is.numeric(mz_range) || length(mz_range) != 2L || diff(mz_range) <= 0) {
    stop("mz_range must be two increasing numbers.")
  }
  style <- match.arg(style)
  points <- get_plot_region_from_file(
    infile, scan_range[[1L]], scan_range[[2L]], mz_range[[1L]], mz_range[[2L]],
    as.integer(ms_level)
  )
  if (length(points) == 0L) stop("No points were found in the requested region.")
  outfile <- .asari_public_output_file(output)
  grDevices::pdf(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  if (style == "scatter") plot_scatter_map_region(points)
  else if (style == "double") double_scatter_map_region(points)
  else with_line_scatter_map_region(points)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  normalizePath(outfile, mustWork = TRUE)
}

#' Perform PCA on feature table sample intensity and plot
#'
#' @param feature_table an explicit asari feature table TSV.
#' @param output Optional PDF output path; `NULL` uses the current drawing device.
#' @return PC1 and PC2 coordinate matrix of the sample.
#' @export
asari_pca <- function(feature_table, output = NULL) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  opened <- !is.null(output)
  if (opened) {
    outfile <- .asari_public_output_file(output)
    grDevices::pdf(outfile)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  scores <- pca_ftable(path)
  if (opened) {
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
    attr(scores, "output_file") <- normalizePath(outfile, mustWork = TRUE)
  }
  scores
}

#' Plot sample summary intensity versus peak_area
#'
#' @param feature_table an explicit asari feature table TSV.
#' @param output explicit PDF output path.
#' @return PDF absolute path.
#' @export
asari_plot_correlations <- function(feature_table, output) {
  path <- path.expand(feature_table)
  if (!file.exists(path) || dir.exists(path)) stop("Feature table does not exist: ", path)
  outfile <- .asari_public_output_file(output)
  grDevices::pdf(outfile)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_correlations(path)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  normalizePath(outfile, mustWork = TRUE)
}
