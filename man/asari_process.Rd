\name{asari_process}
\alias{asari_process}
\title{Process centroided mzML data with asariR}
\usage{
asari_process(
  input,
  output = "asari-output",
  project_name = "asari_project",
  mode = NULL,
  ppm = NULL,
  multicores = NULL,
  rt_align = NULL,
  database_mode = NULL,
  recursive = FALSE,
  parameters = list()
)
}
\arguments{
  \item{input}{One or more explicit \code{.mzML} files, or directories containing them.}
  \item{output}{Parent directory in which a timestamped project is created.}
  \item{project_name}{Short project name used in the output directory name.}
  \item{mode}{Optional ionization mode: \code{"pos"} or \code{"neg"}.}
  \item{ppm}{Optional positive m/z tolerance in parts per million.}
  \item{multicores}{Optional positive number of worker processes.}
  \item{rt_align}{Optional logical value controlling retention-time alignment.}
  \item{database_mode}{Storage mode. The convenience interface defaults to \code{"memory"}.}
  \item{recursive}{Whether directory input should be searched recursively.}
  \item{parameters}{Named list of advanced asari parameter overrides.}
}
\value{
An \code{asari_result} containing the project directory, feature-table paths,
input files, and effective parameters.
}
\description{
This is the main user-facing interface. It validates explicit input paths,
runs the translated asari workflow, and reports the main output tables.
}
\examples{
\dontrun{
result <- asari_process(
  input = "/absolute/path/to/mzML-files",
  output = "/absolute/path/to/results",
  multicores = 1
)
}
}
