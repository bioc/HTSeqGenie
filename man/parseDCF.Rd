\name{parseDCF}
\alias{parseDCF}
\title{Read and parse a configuration file}
\usage{
parseDCF(filename)
}
\arguments{
  \item{filename}{File name}
}
\value{
Named list
}
\description{
From a file like x1: y1 x2: y2 extract field, using the
rules: - split on ':' - first element of split id name of
parameter, second is value - trailing whitespaces (tabs and
spaces) are removed - comments (text flow starting with #)
are removed
}
\keyword{internal}

