% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countGenomicFeatures.R
\name{countFeatures}
\alias{countFeatures}
\title{Count RNA-Seq Pipeline Genomic Features}
\usage{
countFeatures(reads, features)
}
\arguments{
\item{reads}{GRangesList object of interval, usually where reads aligned}

\item{features}{A list of genome annotations as GRangesList}
}
\value{
A list of counts by feature
}
\description{
Given GRanges, counts number of hits by gene, exon, intergenic, etc
}
\details{
Given a GRanges object, this function performs an overlap against
a previously created set of genomic regions. These genomic
regions include genes, coding portions of genes (CDS),
exons, intergenic regions, and exon groups (which
contain two or more exons)
}
\author{
Cory Barr
}
\keyword{internal}

