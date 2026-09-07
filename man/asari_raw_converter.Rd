\name{asari_raw_converter}
\alias{asari_raw_converter}
\alias{asari_install_raw_converter}
\alias{asari_convert_raw}
\alias{asari_uninstall_raw_converter}
\title{Install and use ThermoRawFileParser}
\usage{
asari_install_raw_converter(version = "1.4.5")

asari_convert_raw(
  input,
  output_dir = NULL,
  version = "1.4.5",
  multicores = 1L,
  dry_run = FALSE
)

asari_uninstall_raw_converter(version = "1.4.5")
}
\arguments{
\item{version}{ThermoRawFileParser version.}
\item{input}{One or more explicit Thermo RAW file paths.}
\item{output_dir}{Output directory, or \code{NULL} for each input directory.}
\item{multicores}{Number of conversion workers.}
\item{dry_run}{Whether to return commands without running conversion. The converter must already be installed.}
}
\value{
The installation directory, a conversion-command result list, or invisible
\code{NULL} for uninstallation.
}
\description{
Public wrappers around the translated ThermoRawFileParser installer and batch
converter. Installing the converter requires internet access. macOS and Linux
normally require a compatible Mono runtime.
}
