% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coabundance.R
\name{coabundance}
\alias{coabundance}
\title{Create a coabundance object}
\usage{
coabundance(
  cor_res,
  edges,
  nodes = NULL,
  method = NULL,
  max_pval = NULL,
  min_abs_estimate = NULL,
  ...
)
}
\arguments{
\item{cor_res}{result object originally describing the relationships (e.g. object of class \code{rcorr})}

\item{edges}{tibble with columns from and to describing edges}

\item{method}{character of correlation method used}
}
\description{
Create a coabundance object
}
