\name{sclapply}
\alias{sclapply}
\title{Scheduled parallel processing}
\usage{
sclapply(inext, fun, max.parallel.jobs, ..., stop.onfail = TRUE,
  tracefun = NULL, tracefun.period = 60)
}
\arguments{
  \item{inext}{A function (without argument) returning an
  object to process; NULL if none left; this function is
  run in the main thread}

  \item{fun}{Function to process the object returned by
  inext; this function is run in children thread}

  \item{max.parallel.jobs}{Number of jobs to start in
  parallel}

  \item{...}{Further arguments passed to fun}

  \item{stop.onfail}{Throw error if one}

  \item{tracefun}{Callback function that will be executed
  in a separate thread}

  \item{tracefun.period}{Time intervall between calls to
  tracefun}
}
\value{
Return value of applied function
}
\description{
Scheduled parallel processing
}
\keyword{internal}

