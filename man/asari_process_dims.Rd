\name{asari_process_dims}
\alias{asari_process_dims}
\title{Process direct-infusion mass spectrometry data}
\usage{
asari_process_dims(
  input,
  output = "asari-output",
  project_name = "asari_dims",
  mode = NULL,
  ppm = NULL,
  multicores = NULL,
  database_mode = NULL,
  recursive = FALSE,
  parameters = list()
)
}
\arguments{
  \item{input}{One or more explicit \code{.mzML} files or directories.}
  \item{output}{Parent directory for the timestamped project.}
  \item{project_name}{Short project name.}
  \item{mode}{Optional ionization mode: \code{"pos"} or \code{"neg"}.}
  \item{ppm}{Optional positive m/z tolerance in parts per million.}
  \item{multicores}{Optional positive number of worker processes.}
  \item{database_mode}{Storage mode: \code{"memory"}, \code{"ondisk"}, or \code{"auto"}.}
  \item{recursive}{Whether to search input directories recursively.}
  \item{parameters}{Named list of advanced parameter overrides.}
}
\value{An \code{asari_result}; the main DIMS result is its \code{full_feature_table}.}
\description{
Runs the translated DIMS workflow with retention-time alignment disabled.
}
