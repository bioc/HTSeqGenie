% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{getObjectFilename}
\alias{getObjectFilename}
\title{Get a filename given a directory and the object name}
\usage{
getObjectFilename(dir_path, object_name)
}
\arguments{
\item{dir_path}{A character string containing a dir path}

\item{object_name}{A character string containing the regular expression matching a filename in dir_path}
}
\value{
A character vector containing an existing filename, stops if 0 or more than 1
}
\description{
Get a filename given a directory and the object name
}
\author{
Gregoire Pau
}
\keyword{internal}

