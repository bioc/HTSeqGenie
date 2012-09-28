\name{parseProgressLog}
\alias{parseProgressLog}
\title{Parse a progress.log file}
\usage{
  parseProgressLog(filename, plot = FALSE,
    main = "Processing times", maxt = 30)
}
\arguments{
  \item{filename}{A character string, pointing to a
  "progres.log" file}

  \item{plot}{A logical, indicating whether to draw a
  barplot. Default is FALSE.}

  \item{main}{A character string, indicating the main label
  of the barplot.}

  \item{maxt}{A numeric, indicating the maximal y-range of
  the barplot. Default is 30.}
}
\value{
  A vector containing processing times (in h) and peakmem
  (in GBi)
}
\description{
  Parse a progress.log file
}
\author{
  Gregoire Pau
}
\keyword{internal}

