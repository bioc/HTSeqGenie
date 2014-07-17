\name{markDuplicates}
\alias{markDuplicates}
\title{markDuplicates}
\usage{
markDuplicates(bamfile, outfile = NULL, path = getOption("picard.path"))
}
\arguments{
  \item{bamfile}{Name of input bam file}

  \item{outfile}{Name of output bam file}

  \item{path}{Full path to MarkDuplicates jar}
}
\value{
Path to output bam file
}
\description{
Mark duplicates in bam
}
\details{
Use MarkDuplicates from PicardTools to mark duplicate
alignments in bam file.
}
\author{
Jens Reeder
}

