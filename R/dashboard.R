# Corresponds to Python asari/dashboard.py: project reading, browsing data preparation and injectable Dashboard services.

# This module defines its own NULL fallback value to avoid relying on the loading order of gcms.R.
.dashboard_or <- function(value, fallback) if (is.null(value)) fallback else value

# Corresponds to epd_convert: split into two dictionaries with feature id and empCpd id as keys.
epd_convert <- function(epd_dict) {
  peaks <- list(); compounds <- list()
  for (key in names(epd_dict)) {
    compound <- epd_dict[[key]]
    ids <- character()
    for (peak in compound$MS1_pseudo_Spectra) {
      id <- as.character(peak$id_number)
      peaks[[id]] <- peak
      ids <- c(ids, id)
    }
    compound$MS1_pseudo_Spectra <- ids
    compounds[[key]] <- compound
  }
  list(peaks, compounds)
}

# Corresponds to read_project: read project, cmap, epd and full/preferred feature tables.
read_project <- function(datadir, load_sample_limit = 20L) {
  datadir <- normalizePath(datadir)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Reading project JSON requires jsonlite.")
  project_desc <- jsonlite::fromJSON(file.path(datadir, "project.json"), simplifyVector = FALSE)
  if (!exists(".samples_read_pickle", mode = "function")) stop("Reading dashboard pickle files requires samples.R.")
  cmap <- .samples_read_pickle(file.path(datadir, "export", "cmap.pickle"))
  epd <- .samples_read_pickle(file.path(datadir, "export", "epd.pickle"))
  read_feature_table <- function(path) {
    table <- utils::read.delim(path, row.names = 1L, check.names = FALSE)
    number <- as.integer(.dashboard_or(project_desc$number_of_samples, 0L))
    if (number > load_sample_limit) table <- table[, seq_len(min(ncol(table), 9L + load_sample_limit)), drop = FALSE]
    table
  }
  list(
    project_desc, cmap, epd,
    read_feature_table(file.path(datadir, "export", "full_Feature_table.tsv")),
    read_feature_table(file.path(datadir, "preferred_Feature_table.tsv"))
  )
}

# Corresponds to plot_xic: returns a scatter description that can be drawn by the front end.
plot_xic <- function(xics, mz_dict, track_id) {
  id <- as.character(track_id)
  list(
    type = "scatter", x = xics$rt, y = xics[[id]],
    title = sprintf("Mass track viewer - m/z %4.4f", mz_dict[[id]]),
    ylabel = "Normalized intensity"
  )
}

# Read list or dict by Python track id.
.dashboard_track <- function(cmap, track_id_number) {
  tracks <- cmap$list_mass_tracks
  key <- as.character(as.integer(track_id_number))
  if (!is.null(names(tracks)) && key %in% names(tracks)) tracks[[key]] else tracks[[as.integer(track_id_number) + 1L]]
}

# Corresponds to cmapplot_mass_tracks: Returns the composite mass track drawing description.
cmapplot_mass_tracks <- function(cmap, rt_list, color, track_id_number) {
  id <- as.integer(track_id_number)
  track <- .dashboard_track(cmap, id)
  list(
    type = "scatter", data = data.frame(rt = rt_list, intensity = unlist(track$intensity)),
    color = color,
    title = sprintf("Mass track viewer - (track_id %d, m/z %4.4f)", id, track$mz),
    ylabel = "Intensity (composite)"
  )
}

# Corresponds to convert_dict_html.
convert_dict_html <- function(d, title = "") {
  info <- list(
    c("id_number: ", d$id_number, " - ", "parent_masstrack_id: ", d$parent_masstrack_id, " - ", "parent_epd_id: ", .dashboard_or(d$parent_epd_id, "")),
    c("snr: ", d$snr, " - ", "peak shape: ", round(d$goodness_fitting, 2), " - ", "cSelectivity: ", round(d$cSelectivity, 2)),
    c("height: ", d$height, " - ", "peak_area: ", d$peak_area),
    c("mz: ", round(d$mz, 4), " - ", "rtime: ", paste0("(", paste(round(c(d$rtime_left_base, d$rtime, d$rtime_right_base), 2), collapse = ", "), ")"))
  )
  paste0(title, paste(vapply(info, function(row) paste0("<ul>", paste(row, collapse = " "), "</ul>"), ""), collapse = ""))
}

# Corresponds to convert_dict_markdown: exclude apex/left_base/right_base.
convert_dict_markdown <- function(d, title = "") {
  keys <- setdiff(names(d), c("apex", "left_base", "right_base"))
  paste0(title, paste0("- ", keys, ": \t", vapply(d[keys], as.character, ""), "\n", collapse = ""), "\n")
}

# Corresponds to track_to_peaks.
track_to_peaks <- function(peakDict) {
  result <- list()
  for (peak in peakDict) {
    key <- as.character(peak$parent_masstrack_id)
    result[[key]] <- c(result[[key]], as.character(peak$id_number))
  }
  result
}

# Corresponds to find_track_by_mz: give priority to the nearest track within 0.1Da, otherwise the global closest.
find_track_by_mz <- function(cmap, rt_list, mz) {
  invisible(rt_list)
  tracks <- unname(cmap$list_mass_tracks)
  candidates <- lapply(tracks, function(track) c(abs(mz - track$mz), track$id_number))
  close <- candidates[vapply(candidates, `[[`, 0, 1L) < 0.1]
  use <- if (length(close)) close else candidates
  use[[order(vapply(use, `[[`, 0, 1L), vapply(use, `[[`, 0, 2L))[[1L]]]][[2L]]
}

# Corresponds to find_a_good_peak: Return the first peak with peak shape and selectivity greater than 0.9.
find_a_good_peak <- function(peakDict) {
  Filter(function(peak) peak$goodness_fitting > 0.9 && peak$cSelectivity > 0.9, peakDict)[[1L]]
}

# Read a Python number dictionary, supporting character keys or R positions.
.dashboard_dict_value <- function(dictionary, key0) {
  key <- as.character(key0)
  if (!is.null(names(dictionary)) && key %in% names(dictionary)) dictionary[[key]] else dictionary[[as.integer(key0) + 1L]]
}

# Corresponds to prepare_rt_alignment: completes the mapping and calculates the RT deviation of each sample.
prepare_rt_alignment <- function(cmap) {
  sample_maps <- lapply(cmap$rt_records, `[[`, "reverse_rt_cal_dict")
  rt_length <- as.integer(cmap$rt_length)
  columns <- list()
  all_rt <- NULL
  for (ii in seq_along(sample_maps)) {
    mapping <- sample_maps[[ii]]
    for (key0 in 0:(rt_length - 1L)) if (is.null(mapping[[as.character(key0)]])) mapping[[as.character(key0)]] <- key0
    keys <- sort(as.integer(names(mapping)))
    rt <- vapply(keys, function(key) .dashboard_dict_value(cmap$dict_scan_rtime, key), 0)
    deviation <- vapply(keys, function(key) {
      mapped <- mapping[[as.character(key)]]
      .dashboard_dict_value(cmap$dict_scan_rtime, mapped) - .dashboard_dict_value(cmap$dict_scan_rtime, key)
    }, 0)
    if (is.null(all_rt)) all_rt <- rt
    columns[[cmap$rt_records[[ii]]$name]] <- deviation
  }
  result <- as.data.frame(columns, check.names = FALSE)
  rownames(result) <- as.character(all_rt)
  result
}

# Corresponds to get_summary_panel: returns summary indicators and graph data required for front-end rendering.
get_summary_panel <- function(project_desc, peakDict, epdDict, Ftable, cmap) {
  bins <- as.integer(sqrt(nrow(Ftable)))
  list(
    description = sprintf(
      "Project retrieved from %s, %d features and %d empirical compounds.",
      project_desc$outdir, length(peakDict), length(epdDict)
    ),
    num_bins = bins,
    mz = Ftable$mz,
    rt_peak_area_sqrt = data.frame(rtime = Ftable$rtime, peak_area_sqrt = sqrt(Ftable$peak_area)),
    log10_snr = log10(Ftable$snr + 1),
    peak_shape = Ftable$goodness_fitting,
    cSelectivity = Ftable$cSelectivity,
    rt_alignment = prepare_rt_alignment(cmap),
    show_rt_legend = cmap$`_number_of_samples_` <= 20L
  )
}

# The following six functions correspond to the internal def of dashboard() and are explicitly removed to facilitate testing one by one.
.dashboard_feature_info_by_feature_id <- function(feature_number, peakDict) {
  peak <- peakDict[[paste0("F", as.integer(feature_number))]]
  if (is.null(peak)) sprintf("<p>Feature info not found - %d.</p>", feature_number) else convert_dict_html(peak)
}
.dashboard_cmapplot_track_by_feature_id <- function(feature_number, peakDict, cmap, rt_list) {
  id <- paste0("F", as.integer(feature_number))
  peak <- peakDict[[id]]
  if (is.null(peak)) stop("Feature not found: ", id)
  plot <- cmapplot_mass_tracks(cmap, rt_list, "blue", peak$parent_masstrack_id)
  plot$highlight <- c(peak$rtime_left_base, peak$rtime_right_base, max(1.3 * peak$height, 1e6))
  plot
}
.dashboard_get_features_by_mz <- function(mz, cmap, rt_list, peakDict, track2peaks) {
  id <- as.character(find_track_by_mz(cmap, rt_list, mz))
  features <- lapply(.dashboard_or(track2peaks[[id]], character()), function(feature) peakDict[[feature]])
  list(id, features)
}
.dashboard_track_info_by_mz <- function(mz, cmap, rt_list, peakDict, track2peaks) {
  found <- .dashboard_get_features_by_mz(mz, cmap, rt_list, peakDict, track2peaks)
  if (length(found[[2L]])) paste0(vapply(found[[2L]], function(peak) paste0("<p>", convert_dict_html(peak, peak$id_number), "</p>"), ""), collapse = "")
  else paste0("<p>No qualified feature found on this mass track - ", found[[1L]], ".</p>")
}
.dashboard_cmapplot_track_by_mz <- function(mz, cmap, rt_list, peakDict, track2peaks) {
  found <- .dashboard_get_features_by_mz(mz, cmap, rt_list, peakDict, track2peaks)
  plot <- cmapplot_mass_tracks(cmap, rt_list, "green", found[[1L]])
  plot$highlights <- lapply(found[[2L]], function(peak) c(peak$rtime - 1, peak$rtime + 1, max(1.3 * peak$height, 1e6)))
  plot
}
.dashboard_features_by_rt_with_lines <- function(rtime, rtime_tol, Ftable) {
  subset <- Ftable[Ftable$rtime >= rtime - rtime_tol & Ftable$rtime <= rtime + rtime_tol, , drop = FALSE]
  data.frame(x = subset$mz, y0 = 0, y1 = subset$peak_area)
}

# Corresponds to dashboard: build browsing status; if server adapter is provided, interactive service is started.
dashboard <- function(project_desc, cmap, epd, Ftable, sample_limit = 20L) {
  invisible(sample_limit)
  cat("//*Asari dashboard*//   Press Control-C to exit.\n")
  converted <- epd_convert(epd)
  peakDict <- converted[[1L]]; epdDict <- converted[[2L]]
  good <- find_a_good_peak(peakDict)
  track2peaks <- track_to_peaks(peakDict)
  rt_list <- vapply(0:(as.integer(cmap$rt_length) - 1L), function(ii) .dashboard_dict_value(cmap$dict_scan_rtime, ii), 0)
  state <- list(
    project_desc = project_desc, cmap = cmap, epd = epdDict,
    peakDict = peakDict, Ftable = Ftable, rt_list = rt_list,
    summary = get_summary_panel(project_desc, peakDict, epdDict, Ftable, cmap),
    default_feature_number = as.integer(sub("^F", "", good$id_number)),
    default_mz = round(good$mz, 4), track2peaks = track2peaks
  )
  state$feature_info_by_feature_id <- function(number) .dashboard_feature_info_by_feature_id(number, peakDict)
  state$cmapplot_track_by_feature_id <- function(number) .dashboard_cmapplot_track_by_feature_id(number, peakDict, cmap, rt_list)
  state$get_features_by_mz <- function(mz) .dashboard_get_features_by_mz(mz, cmap, rt_list, peakDict, track2peaks)
  state$track_info_by_mz <- function(mz) .dashboard_track_info_by_mz(mz, cmap, rt_list, peakDict, track2peaks)
  state$cmapplot_track_by_mz <- function(mz) .dashboard_cmapplot_track_by_mz(mz, cmap, rt_list, peakDict, track2peaks)
  state$features_by_rt_with_lines <- function(rt, tol) .dashboard_features_by_rt_with_lines(rt, tol, Ftable)
  server <- getOption("asariR.dashboard_server")
  if (is.function(server)) server(state)
  invisible(state)
}
