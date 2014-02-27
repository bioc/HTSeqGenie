\name{getEndNumber}
\alias{getEndNumber}
\title{Get Read End Number}
\usage{
getEndNumber(int)
}
\arguments{
  \item{int}{an int from a SAM flag}
}
\value{
{1,2}
}
\description{
Returns the end number of an end from a paired-end read
}
\details{
Given an integer from the BAM flag field, tells which end
it is in a read
}
\author{
Cory Barr
}
\keyword{internal}

