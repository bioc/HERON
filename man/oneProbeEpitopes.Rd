% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epitope.R
\name{oneProbeEpitopes}
\alias{oneProbeEpitopes}
\title{Indicate which epitopes are just one probe.}
\usage{
oneProbeEpitopes(epitope_ids)
}
\arguments{
\item{epitope_ids}{vector of epitope ids}
}
\value{
vector of logical indicating epitopes that are one probe
}
\description{
Indicate which epitopes are just one probe.
}
\examples{
oneProbeEpitopes(c("A_1_1", "B_1_1","C_1_2"))
}
