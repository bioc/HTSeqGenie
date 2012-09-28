\name{createTmpDir}
\alias{createTmpDir}
\title{Create a random directory with prefix}
\usage{
  createTmpDir(prefix = NULL)
}
\arguments{
  \item{prefix}{A string that will preced the directory
  name}
}
\value{
  Name of temporary directory
}
\description{
  Especially for testing code it is very helpful to have a
  temp directory with a defined prefix, so one knows which
  test produced which directory.
}
\keyword{internal}

