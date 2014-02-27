\name{countGenomicFeaturesChunk}
\alias{countGenomicFeaturesChunk}
\title{Count reads by genomic Feature}
\usage{
countGenomicFeaturesChunk(save_dir)
}
\arguments{
  \item{save_dir}{PAth to a pipeline run's save dir}
}
\value{
Nothing
}
\description{
Count reads by genomic Feature
}
\details{
given a BAM-file output from gsnap (with the MD tag), count
hits to exons, genes, transcripts, etc. and quantify
miRNA/ncRNA contaminatino
}
\author{
Cory Barr
}
\keyword{internal}

