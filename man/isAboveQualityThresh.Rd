\name{isAboveQualityThresh}
\alias{isAboveQualityThresh}
\title{Check for high quality reads}
\usage{
  isAboveQualityThresh(reads, minquality, minfrac,
    minlength = NULL)
}
\arguments{
  \item{reads}{A set of reads as ShortReadQ object}

  \item{minquality}{Minimal quality score}

  \item{minfrac}{Fraction of positions that need to be over
  minquality to be considered a good read.}

  \item{minlength}{The minimum read length. Default is no
  filter on read length.}
}
\value{
  A boolean vector indicating whether read is considered
  high quality.
}
\description{
  Checks whether reads have more than a fraction of minFrac
  nucleotides with a score below minquality.
}
\keyword{internal}

