\name{asari_dashboard}
\alias{asari_dashboard}
\title{Create an offline asariR project dashboard}
\usage{
asari_dashboard(
  project,
  table = c("preferred", "full"),
  output = NULL,
  open = FALSE,
  max_rows = 5000L
)
}
\arguments{
\item{project}{A project directory or an \code{asari_result}.}
\item{table}{The preferred or full feature table.}
\item{output}{Optional explicit HTML output path.}
\item{open}{Whether to open the generated file in the default browser.}
\item{max_rows}{Maximum number of feature rows embedded in the HTML table.}
}
\value{The absolute path of the generated HTML file.}
\description{
Creates a self-contained searchable HTML summary from \code{project.json} and
the exported feature tables. It does not require Shiny or Python pickle files.
}
