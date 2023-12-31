% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{HERONSequenceDataSet-class}
\alias{HERONSequenceDataSet-class}
\alias{.HERONSequenceDataSet}
\alias{HERONSequenceDataSet}
\title{HERONSequenceDataSet object and constructors}
\usage{
HERONSequenceDataSet(exprs, ...)
}
\arguments{
\item{exprs}{binding values with rows as sequences and columns as samples}

\item{...}{arguments provided to \code{SummarizedExperiment}, including
metadata

metadata can contain a probe DataFrame, that maps sequences
(column PROBE_SEQUENCE) to probe identifiers ( column PROBE_ID)}
}
\value{
HERONSequenceDataSet object
}
\description{
\code{HERONSequenceDataSet} is a subclass of \code{SummarizedExperiment},
used to store the expression values, intermediate calculations, and
results of a differential binding code on the seeuqnce-level.
}
\examples{
exprs <- matrix(seq_len(100),ncol=4)
colnames(exprs) <- c("C1", "C2", "C3", "C4")
sds <- HERONSequenceDataSet(exprs = exprs)
}
