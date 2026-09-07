# Corresponds to the mass operation function of asari/mass_functions.py.
# This file places m/z distance checks, selectivity scoring, clustering, and seed-based grouping logic.
# Compatibility convention: Python 0-based position is converted to R 1-based, tuple return value is converted to named list,
# Python None is converted to NA_integer_; Python calculations, ordering, and boundary rules are preserved under valid scientific input.

# Flatten a list of binary tuples and retain each value exactly once.
#
# The order of traversal is the same as in the Python version: all first elements are collected first, then all second elements.
# R's unique() preserves the order of first occurrence, so the return order is more stable than Python's set().
flatten_tuplelist <- function(tuple_list) {
  if (!is.list(tuple_list)) {
    stop("tuple_list must be a list.", call. = FALSE)
  }
  if (length(tuple_list) == 0L) {
    return(NULL)
  }
  if (any(vapply(tuple_list, length, integer(1)) < 2L)) {
    stop("Each tuple must contain at least two elements.", call. = FALSE)
  }

  first_elements <- lapply(tuple_list, `[[`, 1L)
  second_elements <- lapply(tuple_list, `[[`, 2L)
  unique(unlist(c(first_elements, second_elements), use.names = FALSE))
}

# Find adjacent values in m/z that are sorted in ascending order and are closer than the ppm tolerance.
#
# The Python version returns (ii, ii - 1) starting from 0; here the corresponding R 1-based position is returned,
# Enables the caller to directly use it to index an R vector or list.
check_close_mzs <- function(mzs, mz_tolerance_ppm = 5) {
  if (!is.numeric(mzs) || any(!is.finite(mzs))) {
    stop("mzs must be a finite numeric vector.", call. = FALSE)
  }
  if (length(mz_tolerance_ppm) != 1L ||
      !is.finite(mz_tolerance_ppm) ||
      mz_tolerance_ppm < 0) {
    stop(
      "mz_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }
  if (is.unsorted(mzs, strictly = FALSE)) {
    stop("mzs must be sorted in ascending order.", call. = FALSE)
  }
  if (length(mzs) < 2L) {
    return(list())
  }

  close_pairs <- list()
  for (ii in 2:length(mzs)) {
    tolerance <- mzs[[ii]] * mz_tolerance_ppm * 1e-6
    distance <- mzs[[ii]] - mzs[[ii - 1L]]
    if (distance < tolerance) {
      close_pairs[[length(close_pairs) + 1L]] <- c(ii, ii - 1L)
    }
  }

  close_pairs
}

# Corresponds to calculate_selectivity internal __sel__: convert ppm interval into 0 to 1 selectivity.
`__sel__` <- function(x, std_ppm = 5) {
  if (x > 100) 1 else if (x < 0.1) 0 else 1 - exp(-x / std_ppm)
}

# Compute mass selectivity (mSelectivity) for each m/z in the sorted m/z list.
#
# Each score lies between 0 and 1: the closer it is to 1, the more distinguishable that m/z is from nearby m/z.
# Implement the calculation method that retains the Python version: first calculate the adjacent ppm distance, and then add the current m/z
# Multiply the selectivity scores of up to two lower and two upper neighbors.
calculate_selectivity <- function(sorted_mz_list, std_ppm = 5) {
  if (!is.numeric(sorted_mz_list) || any(!is.finite(sorted_mz_list))) {
    stop("sorted_mz_list must be a finite numeric vector.", call. = FALSE)
  }
  if (length(sorted_mz_list) <= 3L) {
    stop("sorted_mz_list must contain more than three m/z values.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm <= 0) {
    stop("std_ppm must be one finite, positive number.", call. = FALSE)
  }

  # The internal def has been explicitly removed, and the current std_ppm is fixed here.
  selectivity_from_ppm <- function(ppm_distance) `__sel__`(ppm_distance, std_ppm)

  # Python does not force verification of ascending order; when the merged mass track is appended to the end after extraction, the landmark subset may be temporarily non-ascending.
  # At this point, the negative ppm interval is converted to 0 via __sel__, which is exactly the behavior that Python's actual process relies on.
  mz_count <- length(sorted_mz_list)
  ppm_distances <- 1e6 *
    diff(sorted_mz_list) /
    sorted_mz_list[-mz_count]

  # The first two m/z have only incomplete left-hand neighbors and are calculated separately by the Python version.
  selectivities <- c(
    selectivity_from_ppm(ppm_distances[[1L]]) *
      selectivity_from_ppm(ppm_distances[[1L]] + ppm_distances[[2L]]),
    selectivity_from_ppm(ppm_distances[[1L]]) *
      selectivity_from_ppm(ppm_distances[[2L]]) *
      selectivity_from_ppm(ppm_distances[[2L]] + ppm_distances[[3L]])
  )

  # The middle m/z considers both the two lower and the two upper neighbors.
  if (mz_count > 4L) {
    for (ii in 3:(mz_count - 2L)) {
      selectivities <- c(
        selectivities,
        selectivity_from_ppm(
          ppm_distances[[ii - 2L]] + ppm_distances[[ii - 1L]]
        ) *
          selectivity_from_ppm(ppm_distances[[ii - 1L]]) *
          selectivity_from_ppm(ppm_distances[[ii]]) *
          selectivity_from_ppm(
            ppm_distances[[ii]] + ppm_distances[[ii + 1L]]
          )
      )
    }
  }

  # The last two m/z have only incomplete right-hand neighbors and are calculated separately by the Python version.
  distance_count <- length(ppm_distances)
  selectivities <- c(
    selectivities,
    selectivity_from_ppm(
      ppm_distances[[distance_count - 2L]] +
        ppm_distances[[distance_count - 1L]]
    ) *
      selectivity_from_ppm(ppm_distances[[distance_count - 1L]]) *
      selectivity_from_ppm(ppm_distances[[distance_count]]),
    selectivity_from_ppm(
      ppm_distances[[distance_count - 1L]] +
        ppm_distances[[distance_count]]
    ) *
      selectivity_from_ppm(ppm_distances[[distance_count]])
  )

  selectivities
}

# Finds a one-to-one match between two m/z lists without adjacent candidate conflicts.
#
# The function first retains the list source and original position of each m/z, and then merges and sorts the two lists.
# Only if they are from different lists, the distance is less than the ppm tolerance, and the next adjacent value does not cause a conflict.
# Only the current two m/z will be recorded as a match. The R version returns 1-based original position.
mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  list1_count <- length(list1)
  list2_count <- length(list2)
  combined_count <- list1_count + list2_count

  if (combined_count < 2L) {
    return(list(mapped = list(), ratio_deltas = numeric()))
  }

  combined_mz <- c(list1, list2)
  list_origin <- c(rep(1L, list1_count), rep(2L, list2_count))
  original_index <- c(seq_along(list1), seq_along(list2))

  # Python tuples are sorted using m/z, list source, and raw index in order.
  sorted_order <- order(combined_mz, list_origin, original_index)
  combined_mz <- combined_mz[sorted_order]
  list_origin <- list_origin[sorted_order]
  original_index <- original_index[sorted_order]

  # Appends a sentinel value consistent with the Python version so that the last real element is also checked for the next distance.
  combined_mz <- c(combined_mz, 999999)
  list_origin <- c(list_origin, 2L)
  original_index <- c(original_index, NA_integer_)

  mapped <- list()
  ratio_deltas <- numeric()

  for (ii in 2:combined_count) {
    if (list_origin[[ii]] != list_origin[[ii - 1L]]) {
      tolerance <- combined_mz[[ii]] * std_ppm * 1e-6
      distance <- combined_mz[[ii]] - combined_mz[[ii - 1L]]

      if (distance < tolerance &&
          combined_mz[[ii + 1L]] - combined_mz[[ii]] > tolerance) {
        if (list_origin[[ii]] > list_origin[[ii - 1L]]) {
          mapped[[length(mapped) + 1L]] <- c(
            original_index[[ii - 1L]],
            original_index[[ii]]
          )
          ratio_deltas <- c(ratio_deltas, distance / combined_mz[[ii]])
        } else {
          mapped[[length(mapped) + 1L]] <- c(
            original_index[[ii]],
            original_index[[ii - 1L]]
          )
          ratio_deltas <- c(ratio_deltas, -distance / combined_mz[[ii]])
        }
      }
    }
  }

  list(mapped = mapped, ratio_deltas = ratio_deltas)
}

# Complete a one-to-one match between two m/z lists whenever possible.
#
# The function first collects adjacent candidate pairs from different lists that are closer than the ppm tolerance.
# If multiple candidates appear for the same m/z, the pair with the smallest absolute m/z distance is retained.
# Returns the 1-based original position of a successful match, as well as the unmatched positions of each of the two lists.
complete_mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  list1_count <- length(list1)
  list2_count <- length(list2)
  combined_count <- list1_count + list2_count
  candidates <- list()

  if (combined_count >= 2L) {
    combined_mz <- c(list1, list2)
    list_origin <- c(rep(1L, list1_count), rep(2L, list2_count))
    original_index <- c(seq_along(list1), seq_along(list2))

    # Python tuples are sorted using m/z, list source, and raw index in order.
    sorted_order <- order(combined_mz, list_origin, original_index)
    combined_mz <- combined_mz[sorted_order]
    list_origin <- list_origin[sorted_order]
    original_index <- original_index[sorted_order]

    for (ii in 2:combined_count) {
      if (list_origin[[ii]] != list_origin[[ii - 1L]]) {
        tolerance <- combined_mz[[ii]] * std_ppm * 1e-6
        distance <- combined_mz[[ii]] - combined_mz[[ii - 1L]]

        if (distance < tolerance) {
          if (list_origin[[ii]] > list_origin[[ii - 1L]]) {
            candidates[[length(candidates) + 1L]] <- list(
              list1_index = original_index[[ii - 1L]],
              distance = distance,
              list2_index = original_index[[ii]]
            )
          } else {
            candidates[[length(candidates) + 1L]] <- list(
              list1_index = original_index[[ii]],
              distance = distance,
              list2_index = original_index[[ii - 1L]]
            )
          }
        }
      }
    }
  }

  # Python's version of sequential conflict handling: when either original index is shared, the one with the smaller distance wins.
  selected <- list()
  if (length(candidates) > 0L) {
    staged <- candidates[[1L]]

    if (length(candidates) > 1L) {
      for (ii in 2:length(candidates)) {
        current <- candidates[[ii]]
        has_conflict <-
          current$list1_index == staged$list1_index ||
          current$list2_index == staged$list2_index

        if (has_conflict) {
          if (current$distance < staged$distance) {
            staged <- current
          }
        } else {
          selected[[length(selected) + 1L]] <- staged
          staged <- current
        }
      }
    }

    selected[[length(selected) + 1L]] <- staged
  }

  mapped <- lapply(
    selected,
    function(candidate) {
      c(candidate$list1_index, candidate$list2_index)
    }
  )

  mapped_list1 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_list2 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  list(
    mapped = mapped,
    list1_unmapped = setdiff(seq_along(list1), mapped_list1),
    list2_unmapped = setdiff(seq_along(list2), mapped_list2)
  )
}

# Returns all candidates for query m/z within ppm tolerance in the Centurion index.
#
# This is an internal search helper function for all_mass_paired_mapping(). Indexed with
# floor(m/z * 100) is used as a bucket. When querying, check the current bucket and the left and right adjacent buckets.
.find_all_mzmatches_centurion_indexed_list <- function(
    query_mz,
    mz_centurion_tree,
    limit_ppm = 5) {
  query_key <- as.integer(query_mz * 100)
  mz_tolerance <- query_mz * limit_ppm * 1e-6
  results <- list()

  for (bucket_key in (query_key - 1L):(query_key + 1L)) {
    bucket <- mz_centurion_tree[[as.character(bucket_key)]]
    if (is.null(bucket)) {
      next
    }

    for (entry in bucket) {
      if (abs(entry[[1L]] - query_mz) < mz_tolerance) {
        results[[length(results) + 1L]] <- entry
      }
    }
  }

  results
}

# The original Python single-underscore name is retained; the R internal dot-prefixed version is only an early porting alias.
`_find_all_mzmatches_centurion_indexed_list` <- function(
    query_mz, mz_centurion_tree, limit_ppm = 5) {
  .find_all_mzmatches_centurion_indexed_list(
    query_mz, mz_centurion_tree, limit_ppm
  )
}

# Returns all candidate pairs within ppm tolerance between two m/z lists.
#
# Unlike the first two pairing functions, here we do not enforce one-to-one matching, nor do we select only the closest value among multiple candidates.
# A list1 m/z can correspond to multiple list2 m/z, and vice versa. function returns all
# 1-based matching positions, and positions in both lists that are not involved in matching at all.
all_mass_paired_mapping <- function(list1, list2, std_ppm = 5) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }

  # Builds a centurion m/z index of the same concept as mass2chem build_centurion_tree_mzlist().
  mz_centurion_tree <- list()
  for (ii in seq_along(list1)) {
    bucket_key <- as.character(as.integer(list1[[ii]] * 100))
    mz_centurion_tree[[bucket_key]] <- c(
      mz_centurion_tree[[bucket_key]],
      list(c(list1[[ii]], ii))
    )
  }

  mapped <- list()
  for (ii in seq_along(list2)) {
    matches <- .find_all_mzmatches_centurion_indexed_list(
      list2[[ii]],
      mz_centurion_tree,
      limit_ppm = std_ppm
    )

    for (match in matches) {
      mapped[[length(mapped) + 1L]] <- c(as.integer(match[[2L]]), ii)
    }
  }

  mapped_list1 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_list2 <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  list(
    mapped = mapped,
    list1_unmapped = setdiff(seq_along(list1), mapped_list1),
    list2_unmapped = setdiff(seq_along(list2), mapped_list2)
  )
}

# Estimate the overall mass offset of list2 based on reliable m/z matching, correct and re-pair if necessary.
#
# The function first calls mass_paired_mapping() to obtain high-confidence pairings and their relative offsets.
# Then use the average relative offset as the correction ratio of list2. To fully align with the original Python,
# Correction is performed only when the positive offset is strictly greater than correction_tolerance_ppm;
# Negative offsets do not trigger corrections. Returns R 1-based estimated offset ratio before pairing and correction.
mass_paired_mapping_with_correction <- function(
    list1,
    list2,
    std_ppm = 5,
    correction_tolerance_ppm = 1) {
  if (!is.numeric(list1) || any(!is.finite(list1)) ||
      !is.numeric(list2) || any(!is.finite(list2))) {
    stop("list1 and list2 must be finite numeric vectors.", call. = FALSE)
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }
  if (length(correction_tolerance_ppm) != 1L ||
      !is.finite(correction_tolerance_ppm) ||
      correction_tolerance_ppm < 0) {
    stop(
      "correction_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }

  initial_mapping <- mass_paired_mapping(list1, list2, std_ppm)
  mapped <- initial_mapping$mapped
  correction_ratio <- mean(initial_mapping$ratio_deltas)

  if (isTRUE(correction_ratio > correction_tolerance_ppm * 1e-6)) {
    corrected_list2 <- list2 / (1 + correction_ratio)
    mapped <- mass_paired_mapping(
      list1,
      corrected_list2,
      std_ppm
    )$mapped
  }

  list(mapped = mapped, correction_ratio = correction_ratio)
}

# Using the high-confidence landmark m/z as the anchor point, add the m/z list of a new sample to the current reference list.
#
# The function first aligns the landmarks of the reference and new samples, and estimates the overall number when the number of high-confidence pairs is sufficient.
# m/z offset and correct new samples when the absolute value of a positive or negative offset exceeds a threshold. Then align the remaining m/z,
# Update the reference m/z with the average of the paired ends and append the unpaired m/z in the new sample as a new reference row.
# The R version uses 1-based landmarks and mapped locations, using NA_integer_ to represent None in Python.
landmark_guided_mapping <- function(
    REF_reference_mzlist,
    REF_landmarks,
    SM_mzlist,
    SM_landmarks,
    std_ppm = 5,
    correction_tolerance_ppm = 1) {
  if (!is.numeric(REF_reference_mzlist) ||
      any(!is.finite(REF_reference_mzlist)) ||
      !is.numeric(SM_mzlist) ||
      any(!is.finite(SM_mzlist))) {
    stop(
      "REF_reference_mzlist and SM_mzlist must be finite numeric vectors.",
      call. = FALSE
    )
  }
  if (!is.numeric(REF_landmarks) || any(!is.finite(REF_landmarks)) ||
      any(REF_landmarks != as.integer(REF_landmarks)) ||
      any(REF_landmarks < 1L) ||
      any(REF_landmarks > length(REF_reference_mzlist))) {
    stop(
      "REF_landmarks must contain valid 1-based reference positions.",
      call. = FALSE
    )
  }
  if (!is.numeric(SM_landmarks) || any(!is.finite(SM_landmarks)) ||
      any(SM_landmarks != as.integer(SM_landmarks)) ||
      any(SM_landmarks < 1L) ||
      any(SM_landmarks > length(SM_mzlist))) {
    stop(
      "SM_landmarks must contain valid 1-based sample positions.",
      call. = FALSE
    )
  }
  if (length(std_ppm) != 1L || !is.finite(std_ppm) || std_ppm < 0) {
    stop("std_ppm must be one finite, non-negative number.", call. = FALSE)
  }
  if (length(correction_tolerance_ppm) != 1L ||
      !is.finite(correction_tolerance_ppm) ||
      correction_tolerance_ppm < 0) {
    stop(
      "correction_tolerance_ppm must be one finite, non-negative number.",
      call. = FALSE
    )
  }

  REF_landmarks <- as.integer(REF_landmarks)
  SM_landmarks <- as.integer(SM_landmarks)
  original_reference_count <- length(REF_reference_mzlist)
  new_reference_map2 <- rep(NA_integer_, original_reference_count)
  correction_ratio <- NULL

  # Step one: estimate the m/z shift between samples using only high-confidence landmarks.
  reference_anchors <- REF_reference_mzlist[REF_landmarks]
  sample_anchors <- SM_mzlist[SM_landmarks]
  anchor_result <- mass_paired_mapping(
    reference_anchors,
    sample_anchors,
    std_ppm
  )
  anchor_mapping <- anchor_result$mapped

  if (length(anchor_mapping) > 0.2 * original_reference_count) {
    correction_ratio <- mean(anchor_result$ratio_deltas)

    if (abs(correction_ratio) > correction_tolerance_ppm * 1e-6) {
      SM_mzlist <- SM_mzlist / (1 + correction_ratio)
      sample_anchors <- SM_mzlist[SM_landmarks]
      anchor_result <- mass_paired_mapping(
        reference_anchors,
        sample_anchors,
        std_ppm
      )
      anchor_mapping <- anchor_result$mapped
    }
  }

  # Restore positions in the landmark sublist to full m/z list positions.
  mapped <- lapply(
    anchor_mapping,
    function(pair) {
      c(REF_landmarks[[pair[[1L]]]], SM_landmarks[[pair[[2L]]]])
    }
  )
  mapped_reference <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 1L)
  }
  mapped_sample <- if (length(mapped) == 0L) {
    integer()
  } else {
    vapply(mapped, `[[`, integer(1), 2L)
  }

  # Step 2: After excluding landmark matching, match the remaining m/z as much as possible.
  remaining_reference <- setdiff(
    seq_along(REF_reference_mzlist),
    mapped_reference
  )
  remaining_sample <- setdiff(seq_along(SM_mzlist), mapped_sample)
  remaining_result <- complete_mass_paired_mapping(
    REF_reference_mzlist[remaining_reference],
    SM_mzlist[remaining_sample],
    std_ppm
  )
  sample_unmatched <- remaining_sample[remaining_result$list2_unmapped]
  remaining_mapped <- lapply(
    remaining_result$mapped,
    function(pair) {
      c(
        remaining_reference[[pair[[1L]]]],
        remaining_sample[[pair[[2L]]]]
      )
    }
  )
  mapped_pairs <- c(mapped, remaining_mapped)

  message(
    sprintf(
      "    mapped pairs = %d / %d ",
      length(mapped_pairs),
      length(SM_mzlist)
    )
  )

  # On successful pairing, the new sample position is recorded and the reference m/z is updated with the average of both ends of the pair.
  for (pair in mapped_pairs) {
    new_reference_map2[[pair[[1L]]]] <- pair[[2L]]
    REF_reference_mzlist[[pair[[1L]]]] <- 0.5 * (
      REF_reference_mzlist[[pair[[1L]]]] + SM_mzlist[[pair[[2L]]]]
    )
  }

  # m/z that do not have a reference pair in the new sample are appended as new reference rows to the MassGrid.
  for (ii in seq_along(sample_unmatched)) {
    sample_index <- sample_unmatched[[ii]]
    new_reference_map2 <- c(new_reference_map2, sample_index)

    if (sample_index %in% SM_landmarks) {
      REF_landmarks <- c(REF_landmarks, original_reference_count + ii)
    }
  }

  new_reference_mzlist <- c(
    REF_reference_mzlist,
    SM_mzlist[sample_unmatched]
  )

  list(
    new_reference_mzlist = new_reference_mzlist,
    new_reference_map2 = new_reference_map2,
    REF_landmarks = REF_landmarks,
    correction_ratio = correction_ratio
  )
}

# Divides sorted (value, object) into consecutive bins by the current bin median and a user-supplied tolerance function.
#
# Each new value is compared to the median value of the current last bin. The distance is strictly less than
# func_tolerance(value) is put into the current bin, otherwise a new bin is started. Only the object of each tuple is retained when returning.
bin_by_median <- function(List_of_tuples, func_tolerance) {
  if (!is.list(List_of_tuples) || length(List_of_tuples) == 0L) {
    stop("List_of_tuples must be a non-empty list.", call. = FALSE)
  }
  if (any(vapply(List_of_tuples, length, integer(1)) < 2L)) {
    stop("Each tuple must contain a value and an object.", call. = FALSE)
  }
  if (!is.function(func_tolerance)) {
    stop("func_tolerance must be a function.", call. = FALSE)
  }

  bins <- list(list(List_of_tuples[[1L]]))
  if (length(List_of_tuples) > 1L) {
    for (ii in 2:length(List_of_tuples)) {
      current <- List_of_tuples[[ii]]
      current_bin <- bins[[length(bins)]]
      current_median <- stats::median(
        vapply(current_bin, function(tuple) tuple[[1L]], numeric(1))
      )

      if (current[[1L]] - current_median < func_tolerance(current[[1L]])) {
        bins[[length(bins)]][[length(current_bin) + 1L]] <- current
      } else {
        bins[[length(bins) + 1L]] <- list(current)
      }
    }
  }

  lapply(
    bins,
    function(bin) lapply(bin, `[[`, 2L)
  )
}

# Python modules have their own namespace; this alias prevents functions of the same name in tools/merge.R from being obscured in the R flat space.
mass_functions_bin_by_median <- bin_by_median

# Corresponds to gap_divide_mz_cluster internal __divide_by_largest_gap__.
`__divide_by_largest_gap__` <- function(tuple_list) {
  mzs <- vapply(tuple_list, function(tuple) tuple[[1L]], numeric(1))
  split_position <- which.max(diff(mzs)) + 1L
  list(
    tuple_list[seq_len(split_position - 1L)],
    tuple_list[split_position:length(tuple_list)]
  )
}

# Find the largest adjacent gap in the sorted m/z tuple and split the data into two left and right groups from that position.
#
# The mz_tolerance parameter remains in the function interface as in vanilla Python, but the current algorithm does not use it.
# If there are multiple identical maximum intervals, the first one is selected as in numpy.argmax().
gap_divide_mz_cluster <- function(bin_data_tuples, mz_tolerance) {
  if (!is.list(bin_data_tuples) || length(bin_data_tuples) < 2L) {
    stop("bin_data_tuples must contain at least two tuples.", call. = FALSE)
  }
  if (any(vapply(bin_data_tuples, length, integer(1)) < 1L)) {
    stop("Each tuple must contain an m/z value.", call. = FALSE)
  }

  `__divide_by_largest_gap__`(bin_data_tuples)
}

# Replicate the one-dimensional mean filtering of SciPy uniform_filter1d(size, mode="nearest") with nearest boundary extension.
#
# asari passes an integer count to SciPy, SciPy will retain the integer dtype and truncate the decimal part;
# Explicit use of as.integer() here preserves that behavior. Even-numbered windows take one more position to the left.
.uniform_filter1d_nearest_integer <- function(values, size) {
  if (length(values) == 0L || length(size) != 1L || size < 1L) {
    stop("values must be non-empty and size must be positive.", call. = FALSE)
  }

  size <- as.integer(size)
  left_width <- size %/% 2L
  right_width <- size - left_width - 1L
  padded <- c(
    rep(values[[1L]], left_width),
    values,
    rep(values[[length(values)]], right_width)
  )
  cumulative <- c(0, cumsum(padded))
  value_count <- length(values)
  window_sums <- cumulative[(size + 1L):(size + value_count)] -
    cumulative[seq_len(value_count)]

  as.integer(window_sums / size)
}

# Copy scipy.signal.find_peaks(values, distance=...) for peak position selection under the parameters used by asari.
#
# Endpoints are not considered peaks; flat peaks are returned to the midpoint and rounded down. When peak distances conflict, the higher peak is retained first.
# When the heights are equal, the peaks on the right side are processed first in SciPy order, and finally the retained positions are returned in ascending order.
.find_peaks_with_distance <- function(values, distance) {
  if (length(distance) != 1L || !is.finite(distance) || distance < 1) {
    stop("distance must be at least 1.", call. = FALSE)
  }
  value_count <- length(values)
  if (value_count < 3L) {
    return(integer())
  }

  candidates <- integer()
  ii <- 2L
  while (ii <= value_count - 1L) {
    if (values[[ii]] > values[[ii - 1L]]) {
      plateau_end <- ii
      while (plateau_end < value_count &&
             values[[plateau_end + 1L]] == values[[ii]]) {
        plateau_end <- plateau_end + 1L
      }

      if (plateau_end < value_count &&
          values[[plateau_end]] > values[[plateau_end + 1L]]) {
        candidates <- c(candidates, (ii + plateau_end) %/% 2L)
      }
      ii <- plateau_end + 1L
    } else {
      ii <- ii + 1L
    }
  }

  if (length(candidates) <= 1L) {
    return(candidates)
  }

  minimum_distance <- ceiling(distance)
  priority_order <- order(-values[candidates], -candidates)
  kept <- integer()
  for (candidate in candidates[priority_order]) {
    if (length(kept) == 0L ||
        all(abs(candidate - kept) >= minimum_distance)) {
      kept <- c(kept, candidate)
    }
  }

  sort(kept)
}

# Project the m/z tuple to a 1/10000 m/z integer grid and identify mass peak that can serve as clustering seeds.
#
# The function gradually replicates the original Python version: multiply m/z by 10000 and then truncate it to an integer, calculate the frequency of each grid,
# Use nearest mode integer mean filtering, and then use tol4 minimum distance to find peaks. Return values are restored to 0.0001 m/z.
identify_mass_peaks <- function(
    bin_data_tuples,
    mz_tolerance,
    presorted = TRUE) {
  if (!is.list(bin_data_tuples) || length(bin_data_tuples) == 0L) {
    stop("bin_data_tuples must be a non-empty list.", call. = FALSE)
  }
  if (any(vapply(bin_data_tuples, length, integer(1)) < 1L)) {
    stop("Each tuple must contain an m/z value.", call. = FALSE)
  }
  if (length(mz_tolerance) != 1L ||
      !is.finite(mz_tolerance) ||
      mz_tolerance <= 0) {
    stop("mz_tolerance must be one finite, positive number.", call. = FALSE)
  }

  mz4 <- as.integer(
    vapply(bin_data_tuples, function(tuple) tuple[[1L]], numeric(1)) * 10000
  )
  if (!presorted) {
    mz4 <- sort(mz4)
  }

  tol4 <- as.integer(mz_tolerance * 10000)
  if (tol4 < 1L) {
    stop("mz_tolerance must correspond to at least 0.0001 m/z.", call. = FALSE)
  }
  filter_size <- max(2L, as.integer(0.5 * tol4))
  positioned <- seq.int(mz4[[1L]], mz4[[length(mz4)]])
  count_indices <- mz4 - positioned[[1L]] + 1L
  counts <- tabulate(count_indices, nbins = length(positioned))
  filtered_counts <- .uniform_filter1d_nearest_integer(
    counts,
    size = filter_size
  )
  peak_indices <- .find_peaks_with_distance(
    filtered_counts,
    distance = tol4
  )

  0.0001 * positioned[peak_indices]
}

# Assign each tuple to the nearest seed based on the m/z seeds identified by identify_mass_peaks().
#
# When distances are tied, select the earlier seed and preserve the original tuple order within each cluster.
# If the seed is not recognized, use gap_divide_mz_cluster() instead to divide into two groups at the maximum m/z separation as in the Python original.
nn_cluster_by_mz_seeds <- function(
    bin_data_tuples,
    mz_tolerance,
    presorted = TRUE) {
  mz_seeds <- identify_mass_peaks(
    bin_data_tuples,
    mz_tolerance,
    presorted = presorted
  )

  if (length(mz_seeds) > 0L) {
    clusters <- lapply(mz_seeds, function(seed) list())
    for (tuple in bin_data_tuples) {
      distances <- abs(tuple[[1L]] - mz_seeds)
      cluster_index <- which.min(distances)
      clusters[[cluster_index]][[length(clusters[[cluster_index]]) + 1L]] <-
        tuple
    }
    clusters
  } else {
    gap_divide_mz_cluster(bin_data_tuples, mz_tolerance)
  }
}
