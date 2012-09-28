\name{safeUnlink}
\alias{safeUnlink}
\title{safeUnlink}
\usage{
  safeUnlink(path)
}
\arguments{
  \item{path}{A character string indicating which
  file/directory to delete.}
}
\value{
  Nothing
}
\description{
  Symlink-safe file/directory delete function
}
\details{
  Unlike unlink(), safeUnlink() does not follow symlink
  directories for deletion.
}
\author{
  Gregoire Pau
}
\keyword{internal}

