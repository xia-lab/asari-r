# Corresponds to Python asari/tools/extract_targets.py: Export features by target m/z.

# Corresponds to export_targets. The original Python function reference undefined self and filtered_FeatureTable;
# The R version uses registry to explicitly pass in data.frame, and retains the parameters, sample_id, and outfile interfaces.
export_targets <- function(parameters, registry, sample_id, outfile) {
  invisible(sample_id)
  targets <- parameters$target
  if (is.null(targets) || length(targets) == 0L) return(invisible(NULL))
  feature_table <- if (is.data.frame(registry)) registry else registry$filtered_FeatureTable
  if (!is.data.frame(feature_table) || !("mz" %in% names(feature_table))) {
    stop("registry must provide filtered_FeatureTable with an mz column.", call. = FALSE)
  }
  mapping <- all_mass_paired_mapping(
    as.numeric(feature_table$mz), as.numeric(targets), parameters$mz_tolerance_ppm
  )
  unmatched <- mapping$list2_unmapped
  cat(sprintf(
    "\nIn targeted extraction, %d target mz values are not found in this dataset: \n",
    length(unmatched)
  ))
  cat("    ", paste(targets[unmatched], collapse = ", "), "\n")
  pairs <- mapping$mapped
  matched_targets <- if (length(pairs)) {
    vapply(pairs, function(pair) targets[[pair[[2L]]]], 0)
  } else numeric()
  row_index <- if (length(pairs)) vapply(pairs, `[[`, 0L, 1L) else integer()
  targeted <- feature_table[row_index, , drop = FALSE]
  targeted <- cbind(query_target = matched_targets, targeted, stringsAsFactors = FALSE)
  target_outfile <- if (nzchar(outfile)) outfile else file.path(
    parameters$outdir,
    paste0("targeted_extraction__", parameters$output_feature_table)
  )
  utils::write.table(targeted, target_outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  number_samples <- if (!is.null(registry$number_of_samples)) registry$number_of_samples else max(0L, ncol(feature_table) - 11L)
  cat(sprintf(
    "Targeted extraction Feature table (%d x %d) was written to %s.\n\n",
    nrow(targeted), number_samples, target_outfile
  ))
  invisible(NULL)
}
