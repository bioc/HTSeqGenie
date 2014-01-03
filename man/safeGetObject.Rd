\name{safeGetObject}
\alias{safeGetObject}
\title{Safely load a R data file}
\usage{
  safeGetObject(dir_path, object_name)
}
\arguments{
  \item{dir_path}{Save dir of a pipeline run}

  \item{object_name}{object name, can be a regexp}
}
\value{
  loaded object
}
\description{
  Attempts to load a file given by object_name. Bails out
  if none or more than one files match the object name.
}
\keyword{internal}

