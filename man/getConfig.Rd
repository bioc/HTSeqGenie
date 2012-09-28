\name{getConfig}
\alias{getConfig}
\title{Get a configuration parameter}
\usage{
  getConfig(p, stop.ifempty = FALSE)
}
\arguments{
  \item{p}{Name of parameter}

  \item{stop.ifempty}{throw error if value is not set,
  otherwise returns NULL}
}
\value{
  If parameter is missing, return the config list otherwise
  return the value of the parameter name as a character
  string throws an exception if the parameter is not
  present in the config
}
\description{
  Get a configuration parameter
}
\keyword{internal}

