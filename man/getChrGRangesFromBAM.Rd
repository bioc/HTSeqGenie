\name{getChrGRangesFromBAM}
\alias{getChrGRangesFromBAM}
\title{Extract GRanges of Chromosome Lengths from a BAM File}
\usage{
  getChrGRangesFromBAM(bam_file)
}
\arguments{
  \item{bam_file}{bam file to extract GRanges of
  chromosomes from based on header info}
}
\value{
  GRanges
}
\description{
  Creates a GRanges of chromosome length by parsing a BAM
  file header
}
\author{
  Cory Barr
}
\keyword{internal}

