\name{asari_extract_targets}
\alias{asari_extract_targets}
\title{Extract target m/z values from an asari feature table}
\usage{
asari_extract_targets(feature_table, targets, ppm = 5, output = NULL)
}
\arguments{
  \item{feature_table}{An explicit feature-table TSV path or data frame.}
  \item{targets}{One or more positive target m/z values.}
  \item{ppm}{Positive m/z tolerance in parts per million.}
  \item{output}{Optional explicit TSV output path.}
}
\value{Matched feature rows with query target and ppm error columns.}
\description{Returns every feature-target candidate inside the strict ppm window.}
