\name{loadConfig}
\alias{loadConfig}
\title{Load configuration file}
\usage{
  loadConfig(filename)
}
\arguments{
  \item{filename}{Path to configuration file}
}
\value{
  Nothing. Called for its side effect, which is setting the
  global config variable.
}
\description{
  Loads the indicated configuration file. Creates and
  installs a global variable that should be accesssed only
  via getConfig().
}
\keyword{internal}

