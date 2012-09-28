\name{safeExecute}
\alias{safeExecute}
\title{Execute function in try catch with trace function}
\usage{
  safeExecute(expr, memtracer = TRUE)
}
\arguments{
  \item{expr}{Expression to safely execute}

  \item{memtracer}{A boolean, to enable/disable a periodic
  memory tracer. Default is TRUE.}
}
\value{
  Nothing
}
\description{
  Requires the logger to be set
}
\author{
  Gregoire Pau
}
\keyword{internal}

