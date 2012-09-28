\name{writeSummary}
\alias{writeSummary}
\title{Write HTML summary}
\usage{
  writeSummary(dirs, outdir = "Summary", cutoffs)
}
\arguments{
  \item{dirs}{List of pipeline result dirs}

  \item{outdir}{Path to output directory. Does not create
  dir.}

  \item{cutoffs}{list, cutoffs for each plotting/QA
  function}
}
\value{
  Nothing, but writes file
}
\description{
  Write html Summary for list of runs
}
\author{
  Jens Reeder
}
\keyword{internal}

