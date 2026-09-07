\name{asari_analyze}
\alias{asari_analyze}
\alias{asari_extract_mass_tracks}
\alias{print.asari_analysis}
\title{Analyze one mzML file and extract mass-track summaries}
\usage{
asari_analyze(
  input,
  ppm = 5,
  min_intensity = 1000,
  min_timepoints = 6L,
  min_peak_height = 100000,
  parameters = list()
)

asari_extract_mass_tracks(
  input,
  output = NULL,
  ppm = 5,
  min_intensity = 1000,
  min_timepoints = 6L,
  min_peak_height = 100000,
  parameters = list()
)

\method{print}{asari_analysis}(x, ...)
}
\arguments{
  \item{input}{One explicit \code{.mzML} file.}
  \item{ppm}{Positive m/z tolerance in parts per million.}
  \item{min_intensity}{Minimum intensity of retained data points.}
  \item{min_timepoints}{Minimum number of data points in a mass track.}
  \item{min_peak_height}{Minimum peak height.}
  \item{parameters}{Named list of advanced parameter overrides.}
  \item{output}{Optional explicit TSV output path.}
  \item{x}{An \code{asari_analysis} object.}
  \item{...}{Additional print arguments, currently ignored.}
}
\value{
\code{asari_analyze()} returns a structured analysis list.
\code{asari_extract_mass_tracks()} returns a mass-track summary data frame.
}
\description{
Reads one explicit mzML file, extracts mass tracks, and reports scan and
isotope-landmark statistics without connecting to an annotation database.
}
