\name{asari_annotation}
\alias{asari_annotation}
\alias{asari_load_lc_database}
\alias{asari_annotate_lc}
\alias{asari_annotate_gc}
\title{LC-MS and GC-MS annotation interfaces}
\usage{
asari_load_lc_database(database = "hmdb4", cache = TRUE)

asari_annotate_lc(
  feature_table,
  database = "hmdb4",
  output = NULL,
  mode = c("pos", "neg"),
  ppm = 5,
  rt_tolerance = NULL,
  include_adducts = TRUE
)

asari_annotate_gc(
  feature_table,
  kovats,
  library,
  output_dir,
  project_name = "result",
  denovo = FALSE,
  mirror_plots = FALSE,
  ...
)
}
\arguments{
\item{database}{A data frame, explicit database path, or \code{"hmdb4"}.}
\item{cache}{Whether to cache the converted HMDB database as an RDS file.}
\item{feature_table}{An asari feature data frame or explicit TSV path. GC annotation requires a TSV path.}
\item{output}{Optional explicit LC annotation TSV path.}
\item{mode}{Ionization mode, \code{"pos"} or \code{"neg"}.}
\item{ppm}{Mass tolerance in parts per million.}
\item{rt_tolerance}{Optional retention-time tolerance in the input table's units.}
\item{include_adducts}{Whether to generate the default asari adduct masses for neutral compounds.}
\item{kovats}{Explicit Kovats calibration TSV path.}
\item{library}{Explicit MSP or JSON EI library path.}
\item{output_dir}{GC annotation output directory.}
\item{project_name}{Name used in GC output files.}
\item{denovo}{Whether to construct de novo pseudospectra.}
\item{mirror_plots}{Whether to request optional mirror plots.}
\item{...}{Advanced arguments passed to the GC annotation implementation.}
}
\value{
\code{asari_load_lc_database()} returns a normalized data frame.
\code{asari_annotate_lc()} returns candidate feature-compound matches.
\code{asari_annotate_gc()} returns a named list of output paths.
}
\description{
These functions make the translated annotation modules directly callable.
LC annotation performs local mass and optional retention-time candidate matching.
GC annotation uses explicit Kovats standards and an EI spectral library.
}
\details{
The local LC candidate matcher does not claim to reproduce JMS/Khipu empirical
compound grouping. The first read of the original Python HMDB pickle requires a
working Python interpreter; subsequent calls may use the RDS cache.
}
\examples{
features <- data.frame(id_number = "F1", mz = 100.0002, rtime = 10)
standards <- data.frame(id = "S1", name = "standard", mz = 100, rtime = 10)
asari_annotate_lc(features, standards, mode = "pos", ppm = 5)
}
