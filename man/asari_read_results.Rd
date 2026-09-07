\name{asari_read_results}
\alias{asari_read_results}
\alias{print.asari_project_results}
\title{Read an asariR project without dashboard dependencies}
\usage{
asari_read_results(
  project,
  table = c("both", "preferred", "full"),
  max_samples = NULL
)

\method{print}{asari_project_results}(x, ...)
}
\arguments{
  \item{project}{An existing project directory or an \code{asari_result}.}
  \item{table}{Whether to read the preferred table, full table, or both.}
  \item{max_samples}{Optional maximum number of sample abundance columns.}
  \item{x}{An \code{asari_project_results} object.}
  \item{...}{Additional print arguments, currently ignored.}
}
\value{An \code{asari_project_results} containing project metadata and selected tables.}
\description{
Reads the normal project JSON and feature tables. It does not require the
annotation-only CMAP or empirical-compound pickle files used by the dashboard.
}
