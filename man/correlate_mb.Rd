% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlate.R
\name{correlate_mb}
\alias{correlate_mb}
\title{Coabundance analysis using mb as implemented in SpiecEasi}
\usage{
correlate_mb(
  data,
  pulsar.params = list(thresh = 0.05, subsample.ratio = 0.8, ncores =
    getOption("mc.cores"), rep.num = 20, seed = 1337),
  ...
)
}
\arguments{
\item{data}{integer matrix of abundance count data. One sample per row and one taxon per column}

\item{pulsar.params}{list of options passed to \link[SpiecEasi]{pulsar.params}}

\item{...}{further options passed to \link[SpiecEasi]{spiec.easi}}
}
\description{
This is a wrapper arround function \code{SpiecEasi::spiec.easi} with argument \code{method = "mb"}.
}
\references{
\insertRef{spiec_easi}{coabundance}

\insertRef{mb}{coabundance}
}
