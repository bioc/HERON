% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{getProteinLabel}
\alias{getProteinLabel}
\title{Get Protein Label from Probe}
\usage{
getProteinLabel(probes)
}
\arguments{
\item{probes}{vector of probes (i.e. c("A;1", "A;2"))}
}
\value{
vector of strings indicating the protein associated with the
respective probes
}
\description{
Get Protein Label from Probe
}
\examples{
getProteinLabel("A;1")
getProteinLabel("B;2")
getProteinLabel(c("A;1","B;2"))
}
