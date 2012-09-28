\name{removeChunkDir}
\alias{removeChunkDir}
\title{Remove chunk directories}
\usage{
  removeChunkDir()
}
\value{
  Nothing
}
\description{
  Remove chunk directories
}
\details{
  A pipeline run processes the data in small chunks, which
  are eventually combined into the final result.
  Afterwards, this function can be called to remove the
  temporary results per chunk.
}
\author{
  Jens Reeder
}
\keyword{internal}

