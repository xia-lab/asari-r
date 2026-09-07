# Corresponds to Python asari/tools/feature_graph.py: Convert GC feature table into a co-elution relationship graph.

# Read a row of data.frame and convert it into a scalar named list.
.feature_graph_row <- function(data, row) {
  lapply(data[row, , drop = FALSE], function(value) value[[1L]])
}

# Create a variable R environment object corresponding to Python FeatureGraph.
FeatureGraph <- function(ft_path, graph = NULL) FeatureGraph__init__(ft_path, graph)

# Corresponds to FeatureGraph.__init__.
FeatureGraph__init__ <- function(ft_path, graph = NULL) {
  self <- new.env(parent = emptyenv())
  class(self) <- c("FeatureGraph", "environment")
  self$ft_path <- ft_path
  self$df <- if (!is.null(ft_path) && nzchar(ft_path)) {
    utils::read.delim(ft_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else NULL
  self$graph <- if (!is.null(graph)) graph else FeatureGraph_ft_to_graph(self)
  self$clusters <- NULL
  self$reverse_clusters <- NULL
  if (is.null(self$graph)) stop("Either a graph or a feature table path must be provided.")
  if (is.null(self$df)) stop("Either a graph or a feature table path must be provided.")

  # Install a closure with the same name as a Python instance method.
  self$ft_to_graph <- function() FeatureGraph_ft_to_graph(self)
  self$graph_to_ft <- function() FeatureGraph_graph_to_ft(self)
  self$filter_graph <- function(drt = 0.5) FeatureGraph_filter_graph(self, drt)
  self$find_spectral_clusters <- function() FeatureGraph_find_spectral_clusters(self)
  self$extract_fragmentation_spectrum <- function(find_clusters = FALSE, MIN_PEAKS_EXTRACTION = 3L) {
    FeatureGraph_extract_fragmentation_spectrum(self, find_clusters, MIN_PEAKS_EXTRACTION)
  }
  self$map_annotations <- function(matches, to_extract = c("compound_name", "inchikey", "formula")) {
    FeatureGraph_map_annotations(self, matches, to_extract)
  }
  self
}

# Corresponds to the static method ftgraph_from_ft.
FeatureGraph_ftgraph_from_ft <- function(ft_path) FeatureGraph(ft_path)

# Corresponds to the static method ftgraph_from_graph.
FeatureGraph_ftgraph_from_graph <- function(ft_path, graph) FeatureGraph(ft_path, graph)

# Corresponds to the static method metric.
FeatureGraph_metric <- function(x, y) {
  c(abs(as.numeric(x$mz) - as.numeric(y$mz)), abs(as.numeric(x$rtime) - as.numeric(y$rtime)))
}

# Corresponds to ft_to_graph: Create an undirected complete graph containing all feature pairs.
FeatureGraph_ft_to_graph <- function(self) {
  if (is.null(self$df)) return(NULL)
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  node_ids <- vapply(records, function(record) as.character(record$id_number), "")
  nodes <- stats::setNames(records, node_ids)
  edges <- data.frame(u = character(), v = character(), dmz = numeric(), drt = numeric())
  if (length(records) > 1L) {
    pairs <- utils::combn(seq_along(records), 2L)
    edges <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(ii) {
      left <- pairs[1L, ii]
      right <- pairs[2L, ii]
      metric <- FeatureGraph_metric(records[[left]], records[[right]])
      data.frame(
        u = node_ids[[left]], v = node_ids[[right]],
        dmz = metric[[1L]], drt = metric[[2L]], stringsAsFactors = FALSE
      )
    }))
  }
  list(nodes = nodes, edges = edges)
}

# Corresponds to graph_to_ft: restore graph node attributes to data.frame.
FeatureGraph_graph_to_ft <- function(self) {
  records <- unname(self$graph$nodes)
  if (length(records) == 0L || all(lengths(records) == 0L)) return(data.frame())
  keys <- unique(unlist(lapply(records, names), use.names = FALSE))
  rows <- lapply(records, function(record) {
    stats::setNames(lapply(keys, function(key) record[[key]] %||% NA), keys)
  })
  as.data.frame(do.call(rbind, lapply(rows, as.data.frame)), stringsAsFactors = FALSE)
}

# Null substitution operator private to this file.
`%||%` <- function(value, fallback) if (is.null(value)) fallback else value

# Corresponds to filter_graph: only retain the edges and their endpoints where drt is strictly less than the threshold.
FeatureGraph_filter_graph <- function(self, drt = 0.5) {
  selected <- self$graph$edges[self$graph$edges$drt < drt, , drop = FALSE]
  node_ids <- unique(c(selected$u, selected$v))
  # NetworkX add_edges_from only takes (u, v), so the filtered graph does not retain the original edge attributes.
  filtered_edges <- if (nrow(selected)) {
    data.frame(u = selected$u, v = selected$v, dmz = NA_real_, drt = NA_real_)
  } else data.frame(u = character(), v = character(), dmz = numeric(), drt = numeric())
  filtered_nodes <- stats::setNames(rep(list(list()), length(node_ids)), node_ids)
  FeatureGraph_ftgraph_from_graph(
    self$ft_path,
    list(nodes = filtered_nodes, edges = filtered_edges)
  )
}

# Calculate the connected components of the undirected graph and return a list of node id character vectors.
.feature_graph_components <- function(graph) {
  nodes <- names(graph$nodes)
  components <- list()
  unseen <- nodes
  while (length(unseen) > 0L) {
    queue <- unseen[[1L]]
    component <- character()
    while (length(queue) > 0L) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      if (current %in% component) next
      component <- c(component, current)
      edges <- graph$edges
      neighbors <- unique(c(edges$v[edges$u == current], edges$u[edges$v == current]))
      queue <- unique(c(queue, setdiff(neighbors, component)))
    }
    components[[length(components) + 1L]] <- component
    unseen <- setdiff(unseen, component)
  }
  components
}

# Corresponds to find_spectral_clusters: the connection component number starts from 0.
FeatureGraph_find_spectral_clusters <- function(self) {
  if (is.null(self$graph)) stop("Must have a graph to find spectral clusters.")
  components <- .feature_graph_components(self$graph)
  cluster_to_feature <- list()
  feature_to_cluster <- list()
  for (ii in seq_along(components)) {
    cluster_id <- ii - 1L
    cluster_to_feature[[as.character(cluster_id)]] <- components[[ii]]
    for (feature in components[[ii]]) {
      self$graph$nodes[[feature]]$clique_id <- cluster_id
      feature_to_cluster[[feature]] <- cluster_id
    }
  }
  cat(sprintf("Clustering finds: %d spectral clusters!\n", length(cluster_to_feature)))
  self$clusters <- cluster_to_feature
  self$reverse_clusters <- feature_to_cluster
  invisible(NULL)
}

# Create lightweight Spectrum-compatible objects.
.feature_graph_spectrum <- function(mzs, intensities, metadata) {
  list(
    mz = as.numeric(mzs), intensities = as.numeric(intensities),
    metadata = metadata,
    peaks = cbind(mz = as.numeric(mzs), intensity = as.numeric(intensities))
  )
}

# Default spectrum filtering/normalization; users can inject matchms equivalent backend through options.
.feature_graph_process_spectrum <- function(spectrum) {
  filter <- getOption("asariR.spectrum_default_filters")
  normalize <- getOption("asariR.spectrum_normalize")
  if (is.function(filter)) spectrum <- filter(spectrum)
  if (is.function(normalize)) return(normalize(spectrum))
  maximum <- max(spectrum$intensities)
  if (is.finite(maximum) && maximum > 0) {
    spectrum$intensities <- spectrum$intensities / maximum
    spectrum$peaks[, 2L] <- spectrum$intensities
  }
  spectrum
}

# Corresponds to extract_fragmentation_spectrum: Each cluster selects the sample with the highest total intensity and sufficient number of peaks.
FeatureGraph_extract_fragmentation_spectrum <- function(
    self, find_clusters = FALSE, MIN_PEAKS_EXTRACTION = 3L) {
  if (isTRUE(find_clusters)) self$find_spectral_clusters()
  else if (is.null(self$clusters)) stop("Must find clusters before extracting fragmentation spectra.")
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  feature_table <- stats::setNames(records, vapply(records, function(x) as.character(x$id_number), ""))
  samples <- if (ncol(self$df) > 11L) names(self$df)[12:ncol(self$df)] else character()
  processed <- list()
  for (cluster_id in names(self$clusters)) {
    feature_ids <- self$clusters[[cluster_id]]
    spectra <- list()
    for (sample in samples) {
      values <- lapply(feature_ids, function(feature) {
        c(as.numeric(feature_table[[feature]]$mz), as.numeric(feature_table[[feature]][[sample]]))
      })
      values <- values[order(vapply(values, `[[`, 0, 1L), vapply(values, `[[`, 0, 2L))]
      mzs <- vapply(values, `[[`, 0, 1L)
      intensities <- vapply(values, `[[`, 0, 2L)
      if (sum(intensities) > 0) {
        spectra[[length(spectra) + 1L]] <- .feature_graph_spectrum(
          mzs, intensities, list(cluster_id = as.integer(cluster_id), sample = sample)
        )
      }
    }
    if (length(spectra) == 1L) {
      # The original version of Python only assigns values but does not append, retaining this actual behavior.
      selected_spectrum <- spectra[[1L]]
      invisible(selected_spectrum)
    } else if (length(spectra) > 1L) {
      totals <- vapply(spectra, function(spectrum) sum(spectrum$intensities), 0)
      for (index in order(-totals, seq_along(totals))) {
        selected <- .feature_graph_process_spectrum(spectra[[index]])
        if (nrow(selected$peaks) >= MIN_PEAKS_EXTRACTION) {
          processed[[length(processed) + 1L]] <- selected
          break
        }
      }
    }
  }
  processed
}

# Corresponds to map_annotations: propagate the cluster matching to each feature and output the annotation table.
FeatureGraph_map_annotations <- function(
    self, matches, to_extract = c("compound_name", "inchikey", "formula")) {
  records <- lapply(seq_len(nrow(self$df)), function(ii) .feature_graph_row(self$df, ii))
  feature_table <- stats::setNames(records, vapply(records, function(x) as.character(x$id_number), ""))
  for (match in matches) {
    cluster_id <- as.character(match$extract$metadata$cluster_id)
    annotation <- unlist(lapply(to_extract, function(field) {
      c(match$library$metadata[[field]] %||% "Not Found", match$similarity, match$match_peaks)
    }), use.names = FALSE)
    for (feature in self$clusters[[cluster_id]]) {
      existing <- feature_table[[feature]]$annotations %||% list()
      # Python indentation causes the same final annotation to be appended once per field.
      feature_table[[feature]]$annotations <- c(
        existing, rep(list(annotation), length(to_extract))
      )
    }
  }
  # The original column maintains the data.frame type, and the annotation is used as a list-column alone to prevent rbind from misconverting the embedded list into characters.
  ids <- as.character(self$df$id_number)
  self$df$annotations <- I(lapply(ids, function(id) {
    feature_table[[id]]$annotations %||% list()
  }))
  output <- self$df
  if ("annotations" %in% names(output)) {
    output$annotations <- vapply(output$annotations, function(value) paste(capture.output(dput(value)), collapse = ""), "")
  }
  outfile <- gsub(".tsv", "_annotated_gc_beta.tsv", self$ft_path, fixed = TRUE)
  utils::write.table(output, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(NULL)
}
