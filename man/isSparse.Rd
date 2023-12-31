% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.R
\name{isSparse}
\alias{isSparse}
\title{isSparse}
\usage{
isSparse(cov, threshold = 0.1)
}
\arguments{
\item{cov}{A cov object as SimpleRleList}

\item{threshold}{Fraction of number of runs over total length}
}
\value{
Boolean whether this object is dense or sparse
}
\description{
Check coverage for sparseness
}
\details{
Some Rle related operations become very slow when they are dealing with
data that violates their sparseness assumption. This method provides an estimate about
whether the data is dense or sparse. More precicely it checks
if the fraction of the number of runs over the total length is smaller than a threshold
}
\author{
Jens Reeder
}

