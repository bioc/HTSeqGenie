\name{setUpDirs}
\alias{setUpDirs}
\title{Create output directory and subdirectories for sequencing pipeline analysis outputs}
\usage{
  setUpDirs(save_dir, overwrite = "never")
}
\arguments{
  \item{save_dir}{path to the directory that will contain
  all needed subdirectories}

  \item{overwrite}{A character string: never (default),
  erase, overwrite}
}
\value{
  Nothing. Called for its side effects
}
\description{
  Creates a directory with all needed subdirectories for
  pipeline outputs
}
\author{
  Cory Barr, Jens Reeder
}
\keyword{internal}

