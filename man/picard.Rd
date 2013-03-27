\name{picard}
\alias{picard}
\title{picard}
\usage{
  picard(tool, ..., path = getOption("picard.path"))
}
\arguments{
  \item{tool}{Name of the Picard Tool, e.g. MarkDuplicates}

  \item{...}{Arguments forwarded to the picard tool}

  \item{path}{full path to the picard tool jar file.}
}
\value{
  Nothing
}
\description{
  Generic function to call all picard command line java
  tools
}
\author{
  Jens Reeder, Michael Lawrence
}
\keyword{internal}

