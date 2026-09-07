\name{asari_interfaces}
\alias{asari_interfaces}
\alias{asari_gui}
\alias{asari_cli}
\title{Interactive and command-line interfaces}
\usage{
asari_gui(
  input = NULL,
  output = "asari-output",
  project_name = "asari_project",
  mode = c("pos", "neg"),
  workflow = c("LC", "GC", "DIMS"),
  ppm = 5,
  multicores = 1L,
  run = TRUE,
  parameters = list()
)

asari_cli(argv = commandArgs(trailingOnly = TRUE))
}
\arguments{
\item{input}{An explicit mzML file or directory. Interactive sessions may prompt when omitted.}
\item{output}{Output root directory.}
\item{project_name}{Project name.}
\item{mode}{Ionization mode.}
\item{workflow}{LC, GC, or DIMS workflow.}
\item{ppm}{Mass tolerance in parts per million.}
\item{multicores}{Number of worker processes.}
\item{run}{Whether to run processing or only return the prepared configuration.}
\item{parameters}{Named advanced parameter overrides.}
\item{argv}{Command-line argument vector.}
}
\value{
\code{asari_gui()} returns a configuration or processing result.
\code{asari_cli()} invisibly returns the selected command's result.
}
\description{
Provides a dependency-light interactive R/RStudio wizard and a stable Terminal
entry point for process, analyze, annotate, visualization, and workflow listing.
}
