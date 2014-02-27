\name{filterByLength}
\alias{filterByLength}
\title{Filter reads by length}
\usage{
filterByLength(lreads, minlength = 12, paired = FALSE)
}
\arguments{
  \item{lreads}{A set of reads as ShortReadQ object}

  \item{minlength}{Minimum length}

  \item{paired}{Indicates whether lreads has one of two
  elements}
}
\value{
A boolean vector indicating whether read passes filter
}
\description{
Checks whether reads have at least a length of minlength.
Useful values are zero the rid of empty reads or 12 to
match the gsnap k-mer size.
}
\keyword{internal}

