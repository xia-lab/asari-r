\name{asari_analysis_tools}
\alias{asari_analysis_tools}
\alias{asari_merge_feature_tables}
\alias{asari_feature_graph}
\alias{asari_mirror_plot}
\alias{asari_plot_mz_region}
\alias{asari_pca}
\alias{asari_plot_correlations}
\title{Merge, graph, and plot asariR results}
\usage{
asari_merge_feature_tables(
  feature_tables,
  output,
  ppm = 4,
  rt_tolerance = 15,
  method = c("standard", "legacy"),
  make_feature_id = FALSE
)

asari_feature_graph(feature_table, rt_tolerance = 0.5, output = NULL)

asari_mirror_plot(
  query,
  reference,
  output,
  tolerance = NULL,
  title = "MS Mirror Plot"
)

asari_plot_mz_region(
  input,
  output,
  scan_range,
  mz_range,
  ms_level = 1L,
  style = c("scatter", "double", "line")
)

asari_pca(feature_table, output = NULL)

asari_plot_correlations(feature_table, output)
}
\arguments{
\item{feature_tables}{At least two explicit feature-table TSV paths.}
\item{output}{An explicit TSV or PDF output path, as appropriate.}
\item{ppm}{Mass tolerance in parts per million.}
\item{rt_tolerance}{Retention-time tolerance.}
\item{method}{Standard asari tables or the translated legacy formula-mass format.}
\item{make_feature_id}{Passed to the legacy merge implementation.}
\item{feature_table}{An explicit asari feature-table TSV path.}
\item{query, reference}{Two-column mass and intensity matrices or data frames.}
\item{tolerance}{Optional fragment mass tolerance in Da.}
\item{title}{Plot title.}
\item{input}{One explicit mzML path.}
\item{scan_range}{Two increasing scan numbers defining an open interval.}
\item{mz_range}{Two increasing mass values defining an open interval.}
\item{ms_level}{Mass-spectrometry level.}
\item{style}{Region plot style.}
}
\value{
A merged data frame, a feature-cluster data frame, a PCA score matrix, or an
absolute output-file path, depending on the selected function.
}
\description{
Public wrappers for translated merge, feature-graph, and plotting utilities.
The standard merge groups features by mass and retention time and preserves
source IDs and prefixed sample columns.
}
