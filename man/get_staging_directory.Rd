% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pantry.R
\name{get_staging_directory}
\alias{get_staging_directory}
\title{Get the staging directory for the given dataset}
\usage{
get_staging_directory(
  data_set,
  pantry_config = "~/.pantry_config",
  verbose = TRUE
)
}
\arguments{
\item{data_set}{A tag for the given dataset}

\item{pantry_config}{see get_pantry_config}
}
\value{
'<pantry_config$staging_dir>/<data_set>', creating the directory if necessary
}
\description{
Get the staging directory for the given dataset
}
