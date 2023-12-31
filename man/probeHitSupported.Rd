% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multilevel.R
\name{probeHitSupported}
\alias{probeHitSupported}
\title{Find probe hits with a consecutive probe or another sample}
\usage{
probeHitSupported(hit_mat)
}
\arguments{
\item{hit_mat}{matrix of logical values that indicate a hit with a
TRUE value}
}
\value{
matrix of logical values indicate that the TRUE hit is supported by
a consecutive probe hit in the sample sample or the within another sample
}
\description{
Find probe hits with a consecutive probe or another sample
}
