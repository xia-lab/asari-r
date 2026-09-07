\name{asari_ms2}
\alias{asari_ms2}
\alias{asari_extract_ms2}
\alias{asari_match_ms2}
\alias{asari_search_ms2}
\title{Extract, match, and search MS2 spectra}
\usage{
asari_extract_ms2(
  input,
  output = NULL,
  min_intensity = 1000,
  ms2_peak_limit = 50L
)

asari_match_ms2(
  feature_table,
  ms2_files,
  output,
  rt_tolerance = 30,
  ppm = 5
)

asari_search_ms2(
  spectra,
  library,
  output = NULL,
  ms1_tolerance = 0.01,
  ms2_tolerance = 0.02,
  min_similarity = 0.6,
  method = c("entropy", "cosine"),
  min_intensity = 1000
)
}
\arguments{
\item{input}{One explicit mzML file.}
\item{output}{Optional explicit JSON or TSV output path, as appropriate.}
\item{min_intensity}{Minimum peak intensity used during mzML extraction.}
\item{ms2_peak_limit}{Maximum number of peaks retained per MS2 spectrum, or \code{NULL}.}
\item{feature_table}{Explicit full asari feature-table TSV path.}
\item{ms2_files}{One or more explicit mzML paths.}
\item{rt_tolerance}{Retention-time tolerance in seconds.}
\item{ppm}{Precursor mass tolerance in parts per million for MS1 feature matching.}
\item{spectra}{An mzML path, an \code{asari_extract_ms2()} result, or an MS2 spectrum list.}
\item{library}{An MSP, MGF, or JSON library path, or a parsed library list.}
\item{ms1_tolerance}{Absolute precursor tolerance in Da for library search.}
\item{ms2_tolerance}{Absolute fragment tolerance in Da.}
\item{min_similarity}{Minimum accepted score between zero and one.}
\item{method}{Local entropy-style or cosine scoring.}
}
\value{
Extracted spectra, matched spectra by MS1 feature, or a library-hit data frame.
}
\description{
Public interfaces for the translated MS2 extraction, MS1 matching, and local
spectral-library search code.
}
\details{
When no external entropy backend is configured, \code{method = "entropy"} uses
the package's documented local Jensen-Shannon similarity fallback and is not
guaranteed to be numerically identical to the Python \code{ms_entropy} library.
}
