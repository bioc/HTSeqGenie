\name{storeCoverage}
\alias{storeCoverage}
\title{Save coverage information}
\usage{
  storeCoverage(bam.file, extendReads = FALSE,
    paired_ends = FALSE, fragmentLength = NULL)
}
\arguments{
  \item{bam.file}{Name of bam file.}

  \item{extendReads}{A logical indicating whether reads
  should be extended. If TRUE, reads are extended using
  mate-pair information if paired_ends is TRUE and using
  fragmentLength otherwise. Default is FALSE.}

  \item{paired_ends}{A logical indicating whether reads are
  paired. Default is FALSE.}

  \item{fragmentLength}{An optional numeric, indicating the
  fragment length. If NULL, fragment length is estimated
  from the reads using
  \code{chipseq::estimate.mean.fraglen}. Default is NULL.}
}
\value{
  coverage
}
\description{
  Save coverage information
}
\details{
  Saves the coverage as RData and BigWig
}
\author{
  Cory Barr, Jens Reeder
}
\keyword{internal}

