% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Class_denoiSeq.R
\docType{methods}
\name{setStepSizes}
\alias{setStepSizes}
\alias{setStepSizes,readsData-method}
\title{Generic  for altering the stepSizes slot.}
\usage{
setStepSizes(object, stepSizesval)

\S4method{setStepSizes}{readsData}(object, stepSizesval)
}
\arguments{
\item{object}{a readsData object}

\item{stepSizesval}{A list of step sizes for each of  the  parameters.}
}
\value{
The same readsData object with the stepSizes slot updated.
}
\description{
Updates the value of the stepSizes slot for the readsData object supplied.
}
\section{Methods (by class)}{
\itemize{
\item \code{readsData}: Alters the value of the stepSizes slot of
a readsData object.
}}

\examples{
RD <- new("readsData", counts = ERCC)
ss <- list(N_A = rep(2,  92), N_B = rep(1.5, 92), p = 3e5, f = 3.5e7)
RD <- setStepSizes(RD, ss)
RD@stepSizes
}
