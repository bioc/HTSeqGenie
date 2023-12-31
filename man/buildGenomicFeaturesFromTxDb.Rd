% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildGenomicFeaturesFromTxDb.R
\name{buildGenomicFeaturesFromTxDb}
\alias{buildGenomicFeaturesFromTxDb}
\title{Build genomic features from a TxDb object}
\usage{
buildGenomicFeaturesFromTxDb(txdb)
}
\arguments{
\item{txdb}{A TxDb object.}
}
\value{
A list named list of GRanges objects containing the biological entities to account for.
}
\description{
Build genomic features from a TxDb object
}
\examples{
\dontrun{
  library("TxDb.Hsapiens.UCSC.hg19.knownGene")
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  genomic_features <- buildGenomicFeaturesFromTxDb(txdb)
}
}
\author{
Gregoire Pau
}

