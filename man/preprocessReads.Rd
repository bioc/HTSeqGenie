\name{preprocessReads}
\alias{preprocessReads}
\title{Pipeline preprocessing}
\usage{
preprocessReads()
}
\value{
A named vector containg the path to the preproceesed FastQ
files and a few other statistics
}
\description{
The preprocessing for our NGS pipelines consists of : -
quality filtering - check for adapter contamination -
filtering of rRNA reads - read trimming - shortRead report
generation of surviving reads
}
\details{
These steps are mostly controlled by the global config.
}
\keyword{internal}

