\name{getConfig.integer}
\alias{getConfig.integer}
\title{Check if a config parameter is an integer}
\usage{
getConfig.integer(p, tol = 1e-08, ...)
}
\arguments{
  \item{p}{Name of parameter}

  \item{tol}{Tolerance that controls how far a value can be
  from the next integer.}

  \item{...}{Additinal parameters passed to getConfig()}
}
\value{
Value of parameter as integer
}
\description{
Throws exception if value is no integer
}
\keyword{internal}

