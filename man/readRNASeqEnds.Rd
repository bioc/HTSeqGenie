\name{readRNASeqEnds}
\alias{readRNASeqEnds}
\title{Read Paired End Bam Files}
\usage{
  readRNASeqEnds(bam_file, scan_bam_what = NULL,
    remove.strandness = TRUE)
}
\arguments{
  \item{bam_file}{Path to a bam file}

  \item{scan_bam_what}{additional paramters passed to
  scanBamParam()}

  \item{remove.strandness}{A logical indicating whether
  read strands should be set to "*".}
}
\value{
  GRangesList
}
\description{
  Read paired end BAM files with requested columns from the
  BAM
}
\details{
  This function returns a GRangesList containing
  information about pair-end alignment data.  It is
  essentially Rsamtools::readGAlignmentsFromBam, but it
  uses Rsamtools::scanBam to get information that cannot be
  obtained from readGAlignmentsFromBam.
}
\author{
  Cory Barr
}
\keyword{internal}

