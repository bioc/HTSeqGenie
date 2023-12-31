% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{createTmpDir}
\alias{createTmpDir}
\title{Create a random directory with prefix in R temp dir}
\usage{
createTmpDir(prefix = NULL, dir = tempdir())
}
\arguments{
\item{prefix}{A string that will preceed the directory name}

\item{dir}{Directory where the random dir will be created under.
Defaults to tempdir()}
}
\value{
Name of temporary directory
}
\description{
Especially for testing code it is very helpful to have
a temp directory with a defined prefix, so one knows which test
produced which directory.
}
\keyword{internal}

