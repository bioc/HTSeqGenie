\name{initLog}
\alias{initLog}
\title{Initialize the logger}
\usage{
  initLog(save_dir, debug_level = "INFO")
}
\arguments{
  \item{save_dir}{Save dir of a pipeline run}

  \item{debug_level}{One of INFO, WARN, ERROR, FATAL}
}
\value{
  Log file name
}
\description{
  Setup logging file in {save_dir}/progress.log and log
  sessionInfo and configuration
}
\keyword{internal}

