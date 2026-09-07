\name{asari_quality_tools}
\alias{asari_feature_qc}
\alias{asari_qc_report}
\alias{asari_scan_summary}
\title{Create local feature and mzML quality-control reports}
\usage{
asari_feature_qc(feature_table, output, height = 12, aspect = 0.7)

asari_qc_report(input, output, spikeins = NULL)

asari_scan_summary(
  input,
  output,
  nspec_plot = 50L,
  offset_n = 100L,
  width = 10,
  height = 6,
  title = ""
)
}
\arguments{
  \item{feature_table}{An asari full-feature-table TSV path or data frame.}
  \item{input}{One explicit \code{.mzML} file.}
  \item{output}{An explicit PDF or HTML output path, as appropriate.}
  \item{height, width}{Output plot dimensions.}
  \item{aspect}{Width-to-height ratio of the feature QC PDF.}
  \item{spikeins}{\code{NULL}, a target list, or an explicit JSON file.}
  \item{nspec_plot}{Number of example scans to show.}
  \item{offset_n}{Approximate starting scan offset.}
  \item{title}{Optional scan-summary title.}
}
\value{The absolute path of the generated file.}
\description{
Creates local QC output using mzML scans or an already generated feature table.
No annotation database or spectrum library is contacted.
}
