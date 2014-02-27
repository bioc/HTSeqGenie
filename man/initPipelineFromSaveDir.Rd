\name{initPipelineFromSaveDir}
\alias{initPipelineFromSaveDir}
\title{Init Pipeline environment from previous run}
\usage{
initPipelineFromSaveDir(save_dir, config_update)
}
\arguments{
  \item{save_dir}{Save dir of a previous pipeline run}

  \item{config_update}{List of name value pairs that will
  update the config parameters}
}
\value{
Nothing
}
\description{
Init Pipeline environment from previous run
}
\details{
Loads the config file from a previous run stored in
[save_dir]/logs/config.txt
}
\author{
Gregoire Pau
}
\keyword{internal}

