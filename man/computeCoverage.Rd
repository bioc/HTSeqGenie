\name{computeCoverage}
\alias{computeCoverage}
\title{Compute the coverage vector given a bamfile}
\usage{
  computeCoverage(bamfile, extendReads = FALSE,
    paired_ends = FALSE, fragmentLength = NULL,
    maxFragmentLength = NULL)
}
\arguments{
  \item{bamfile}{A character string indicating the path of
  bam file}

  \item{extendReads}{A logical, indicating whether reads
  should be extended}

  \item{paired_ends}{A logical, indicating whether reads
  are paired}

  \item{fragmentLength}{An integer, indicating the new size
  of reads when extendReads is TRUE and paired_ends is
  FALSE. If NULL, read size is estimated using
  estimate.mean.fraglen from the chipseq package.}

  \item{maxFragmentLength}{An optional integer, specifying
  the maximal size of fragments. Longer fragments will be
  disregarded when computing coverage.}
}
\value{
  A SimpleRleList object containing the coverage
}
\description{
  Compute the coverage vector given a bamfile
}
\author{
  Gregoire Pau
}
\keyword{internal}

