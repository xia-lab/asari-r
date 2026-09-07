#!/usr/bin/env Rscript

# 对同一个明确输入产生的Python和R项目结果做单文件差异诊断。
# 本脚本只读取已经生成的特征表，不运行asari，也不会修改核心算法参数。

usage <- function() {
  paste(
    "Usage:",
    "Rscript scripts/diagnose_feature_parity.R",
    "--python-project /absolute/path/to/python-project",
    "--r-project /absolute/path/to/r-project",
    "--output /absolute/path/to/diagnostic-output",
    "[--mz-ppm 5] [--rt-sec 2]"
  )
}

parse_arguments <- function(arguments) {
  values <- list(mz_ppm = 5, rt_sec = 2)
  index <- 1L
  while (index <= length(arguments)) {
    name <- arguments[[index]]
    if (!startsWith(name, "--") || index == length(arguments)) {
      stop(usage(), call. = FALSE)
    }
    key <- sub("^--", "", name)
    key <- gsub("-", "_", key, fixed = TRUE)
    values[[key]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  required <- c("python_project", "r_project", "output")
  missing <- required[!required %in% names(values)]
  if (length(missing) > 0L) {
    stop("Missing arguments: ", paste(missing, collapse = ", "), "\n", usage(), call. = FALSE)
  }
  values$mz_ppm <- as.numeric(values$mz_ppm)
  values$rt_sec <- as.numeric(values$rt_sec)
  if (!is.finite(values$mz_ppm) || values$mz_ppm <= 0 ||
      !is.finite(values$rt_sec) || values$rt_sec <= 0) {
    stop("mz-ppm and rt-sec must be positive numbers.", call. = FALSE)
  }
  values
}

feature_table_path <- function(project, table = c("preferred", "full")) {
  table <- match.arg(table)
  path <- if (table == "preferred") {
    file.path(project, "preferred_Feature_table.tsv")
  } else {
    file.path(project, "export", "full_Feature_table.tsv")
  }
  if (!file.exists(path)) stop("Feature table does not exist: ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}

read_feature_table <- function(path) {
  table <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c(
    "id_number", "mz", "rtime", "rtime_left_base", "rtime_right_base",
    "peak_area", "cSelectivity", "goodness_fitting", "snr", "detection_counts"
  )
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  table
}

# 输出表已经按Python规则把m/z保留4位、RT保留2位；以这两个公开值定义严格同峰。
feature_key <- function(table) sprintf("%.4f|%.2f", table$mz, table$rtime)

# 对严格同键的峰执行一对一配对；重复键优先选择峰面积最接近的记录。
match_exact_features <- function(left, right) {
  if (nrow(left) == 0L || nrow(right) == 0L) {
    return(data.frame(left_row = integer(), right_row = integer()))
  }
  left_keys <- feature_key(left)
  right_keys <- feature_key(right)
  right_groups <- split(seq_len(nrow(right)), right_keys)
  used <- rep(FALSE, nrow(right))
  pairs <- vector("list", nrow(left))
  pair_count <- 0L

  for (left_row in seq_len(nrow(left))) {
    candidates <- right_groups[[left_keys[[left_row]]]]
    candidates <- candidates[!used[candidates]]
    if (length(candidates) == 0L) next
    area_delta <- abs(log1p(right$peak_area[candidates]) - log1p(left$peak_area[[left_row]]))
    right_row <- candidates[[which.min(area_delta)]]
    used[[right_row]] <- TRUE
    pair_count <- pair_count + 1L
    pairs[[pair_count]] <- c(left_row, right_row)
  }

  if (pair_count == 0L) {
    return(data.frame(left_row = integer(), right_row = integer()))
  }
  matrix_pairs <- do.call(rbind, pairs[seq_len(pair_count)])
  data.frame(
    left_row = as.integer(matrix_pairs[, 1L]),
    right_row = as.integer(matrix_pairs[, 2L])
  )
}

# 把推荐表筛选条件写成可读字符串，直接说明某个完整峰为何没有进入推荐表。
failed_preferred_filters <- function(row) {
  failed <- character()
  if (row$detection_counts[[1L]] <= 0) failed <- c(failed, "detection_counts<=0")
  if (row$snr[[1L]] <= 2) failed <- c(failed, "snr<=2")
  if (row$goodness_fitting[[1L]] <= 0.7) failed <- c(failed, "goodness_fitting<=0.7")
  if (row$cSelectivity[[1L]] <= 0.7) failed <- c(failed, "cSelectivity<=0.7")
  if (length(failed) == 0L) "none" else paste(failed, collapse = ";")
}

# 在另一边完整表中寻找严格同峰；若没有，再按ppm和RT容差寻找最接近峰。
classify_unmatched_features <- function(features, other_full, other_label, mz_ppm, rt_sec) {
  rows <- vector("list", nrow(features))
  other_keys <- feature_key(other_full)

  for (index in seq_len(nrow(features))) {
    feature <- features[index, , drop = FALSE]
    exact <- which(other_keys == feature_key(feature))
    classification <- "not_detected_in_other_full"
    matched_row <- NA_integer_

    if (length(exact) > 0L) {
      area_delta <- abs(log1p(other_full$peak_area[exact]) - log1p(feature$peak_area[[1L]]))
      matched_row <- exact[[which.min(area_delta)]]
      classification <- "exact_in_other_full_but_not_preferred"
    } else {
      mz_tolerance <- feature$mz[[1L]] * mz_ppm * 1e-6
      candidates <- which(
        abs(other_full$mz - feature$mz[[1L]]) <= mz_tolerance &
          abs(other_full$rtime - feature$rtime[[1L]]) <= rt_sec
      )
      if (length(candidates) > 0L) {
        score <- abs(other_full$mz[candidates] - feature$mz[[1L]]) / mz_tolerance +
          abs(other_full$rtime[candidates] - feature$rtime[[1L]]) / rt_sec
        matched_row <- candidates[[which.min(score)]]
        classification <- "near_in_other_full"
      }
    }

    matched <- if (is.na(matched_row)) NULL else other_full[matched_row, , drop = FALSE]
    rows[[index]] <- data.frame(
      source_id = feature$id_number,
      source_mz = feature$mz,
      source_rtime = feature$rtime,
      source_left_base = feature$rtime_left_base,
      source_right_base = feature$rtime_right_base,
      source_peak_area = feature$peak_area,
      source_cSelectivity = feature$cSelectivity,
      source_goodness_fitting = feature$goodness_fitting,
      source_snr = feature$snr,
      source_detection_counts = feature$detection_counts,
      classification = classification,
      other_id = if (is.null(matched)) NA_character_ else as.character(matched$id_number),
      other_mz = if (is.null(matched)) NA_real_ else matched$mz,
      other_rtime = if (is.null(matched)) NA_real_ else matched$rtime,
      other_cSelectivity = if (is.null(matched)) NA_real_ else matched$cSelectivity,
      other_goodness_fitting = if (is.null(matched)) NA_real_ else matched$goodness_fitting,
      other_snr = if (is.null(matched)) NA_real_ else matched$snr,
      other_detection_counts = if (is.null(matched)) NA_real_ else matched$detection_counts,
      mz_delta_ppm = if (is.null(matched)) NA_real_ else
        abs(matched$mz - feature$mz) / feature$mz * 1e6,
      rt_delta_sec = if (is.null(matched)) NA_real_ else abs(matched$rtime - feature$rtime),
      left_base_delta_sec = if (is.null(matched)) NA_real_ else
        matched$rtime_left_base - feature$rtime_left_base,
      right_base_delta_sec = if (is.null(matched)) NA_real_ else
        matched$rtime_right_base - feature$rtime_right_base,
      peak_area_ratio = if (is.null(matched) || feature$peak_area[[1L]] == 0) NA_real_ else
        matched$peak_area / feature$peak_area,
      other_failed_filters = if (is.null(matched)) NA_character_ else
        failed_preferred_filters(matched),
      other_implementation = other_label,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

write_tsv <- function(table, path) {
  utils::write.table(
    table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
  )
}

# 汇总严格同峰记录的连续数值误差，判断一次修改是否真正更接近Python。
summarize_pair_errors <- function(left, right, pairs) {
  if (nrow(pairs) == 0L) return(data.frame())
  left_rows <- left[pairs$left_row, , drop = FALSE]
  right_rows <- right[pairs$right_row, , drop = FALSE]
  definitions <- list(
    peak_area_abs = abs(right_rows$peak_area - left_rows$peak_area),
    peak_area_relative = abs(right_rows$peak_area - left_rows$peak_area) /
      pmax(abs(left_rows$peak_area), 1),
    cSelectivity_abs = abs(right_rows$cSelectivity - left_rows$cSelectivity),
    goodness_fitting_abs = abs(right_rows$goodness_fitting - left_rows$goodness_fitting),
    snr_abs = abs(right_rows$snr - left_rows$snr),
    left_base_abs_sec = abs(right_rows$rtime_left_base - left_rows$rtime_left_base),
    right_base_abs_sec = abs(right_rows$rtime_right_base - left_rows$rtime_right_base)
  )
  do.call(rbind, lapply(names(definitions), function(metric) {
    values <- definitions[[metric]]
    data.frame(
      metric = metric,
      exact_count = sum(values == 0, na.rm = TRUE),
      mean = mean(values, na.rm = TRUE),
      median = stats::median(values, na.rm = TRUE),
      maximum = max(values, na.rm = TRUE)
    )
  }))
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
python_project <- normalizePath(arguments$python_project, mustWork = TRUE)
r_project <- normalizePath(arguments$r_project, mustWork = TRUE)
output <- normalizePath(arguments$output, mustWork = FALSE)
dir.create(output, recursive = TRUE, showWarnings = FALSE)

python_preferred <- read_feature_table(feature_table_path(python_project, "preferred"))
python_full <- read_feature_table(feature_table_path(python_project, "full"))
r_preferred <- read_feature_table(feature_table_path(r_project, "preferred"))
r_full <- read_feature_table(feature_table_path(r_project, "full"))

preferred_pairs <- match_exact_features(python_preferred, r_preferred)
full_pairs <- match_exact_features(python_full, r_full)
python_only_rows <- setdiff(seq_len(nrow(python_preferred)), preferred_pairs$left_row)
r_only_rows <- setdiff(seq_len(nrow(r_preferred)), preferred_pairs$right_row)
python_only <- python_preferred[python_only_rows, , drop = FALSE]
r_only <- r_preferred[r_only_rows, , drop = FALSE]
python_only_full_rows <- setdiff(seq_len(nrow(python_full)), full_pairs$left_row)
r_only_full_rows <- setdiff(seq_len(nrow(r_full)), full_pairs$right_row)
python_only_full <- python_full[python_only_full_rows, , drop = FALSE]
r_only_full <- r_full[r_only_full_rows, , drop = FALSE]

python_diagnosis <- classify_unmatched_features(
  python_only, r_full, "R", arguments$mz_ppm, arguments$rt_sec
)
r_diagnosis <- classify_unmatched_features(
  r_only, python_full, "Python", arguments$mz_ppm, arguments$rt_sec
)
# 完整表差异直接反映峰检测阶段的漏报、多报或近邻偏移，不能只看推荐表数量。
python_full_diagnosis <- classify_unmatched_features(
  python_only_full, r_full, "R", arguments$mz_ppm, arguments$rt_sec
)
r_full_diagnosis <- classify_unmatched_features(
  r_only_full, python_full, "Python", arguments$mz_ppm, arguments$rt_sec
)

classification_counts <- function(table, source) {
  if (nrow(table) == 0L) return(data.frame())
  counts <- as.data.frame(table(table$classification), stringsAsFactors = FALSE)
  names(counts) <- c("classification", "count")
  counts$source <- source
  counts[c("source", "classification", "count")]
}

summary <- data.frame(
  metric = c(
    "python_preferred", "r_preferred", "exact_preferred_pairs",
    "python_only_preferred", "r_only_preferred", "python_full", "r_full",
    "exact_full_pairs", "python_only_full", "r_only_full"
  ),
  value = c(
    nrow(python_preferred), nrow(r_preferred), nrow(preferred_pairs),
    nrow(python_only), nrow(r_only), nrow(python_full), nrow(r_full), nrow(full_pairs),
    nrow(python_only_full), nrow(r_only_full)
  )
)
classifications <- rbind(
  classification_counts(python_diagnosis, "Python_only"),
  classification_counts(r_diagnosis, "R_only")
)

write_tsv(summary, file.path(output, "summary.tsv"))
write_tsv(classifications, file.path(output, "classification_counts.tsv"))
write_tsv(python_diagnosis, file.path(output, "python_only_diagnosis.tsv"))
write_tsv(r_diagnosis, file.path(output, "r_only_diagnosis.tsv"))
write_tsv(
  python_full_diagnosis,
  file.path(output, "python_only_full_diagnosis.tsv")
)
write_tsv(
  r_full_diagnosis,
  file.path(output, "r_only_full_diagnosis.tsv")
)
write_tsv(preferred_pairs, file.path(output, "exact_preferred_pairs.tsv"))
write_tsv(full_pairs, file.path(output, "exact_full_pairs.tsv"))
write_tsv(
  summarize_pair_errors(python_full, r_full, full_pairs),
  file.path(output, "exact_full_metric_errors.tsv")
)

cat("Diagnostic output: ", normalizePath(output, mustWork = TRUE), "\n", sep = "")
cat("Python preferred: ", nrow(python_preferred), "\n", sep = "")
cat("R preferred: ", nrow(r_preferred), "\n", sep = "")
cat("Exact one-to-one preferred pairs: ", nrow(preferred_pairs), "\n", sep = "")
cat("Python-only preferred: ", nrow(python_only), "\n", sep = "")
cat("R-only preferred: ", nrow(r_only), "\n", sep = "")
if (nrow(classifications) > 0L) print(classifications, row.names = FALSE)
