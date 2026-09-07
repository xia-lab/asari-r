\name{asari_compare_features}
\alias{asari_compare_features}
\title{Compare two LC-MS feature tables}
\usage{
asari_compare_features(
  left,
  right,
  ppm = 5,
  rt_tolerance = 5,
  left_rt_unit = c("seconds", "minutes"),
  right_rt_unit = c("seconds", "minutes"),
  verbose = FALSE
)
}
\arguments{
  \item{left, right}{Explicit TSV paths or data frames with \code{mz} and \code{rtime}.}
  \item{ppm}{Positive m/z tolerance in parts per million.}
  \item{rt_tolerance}{Positive retention-time tolerance in seconds.}
  \item{left_rt_unit, right_rt_unit}{Retention-time unit in each input.}
  \item{verbose}{Whether to print the underlying bidirectional matching counts.}
}
\value{A data frame of mutually best matches with m/z and retention-time differences.}
\description{
Chooses the smallest m/z-error candidate in both directions within the ppm
and retention-time windows, and keeps only mutually best pairs.
}
