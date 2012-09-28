\name{runPipelineConfig}
\alias{runPipelineConfig}
\title{Run the NGS analysis pipeline}
\usage{
  runPipelineConfig(config_filename, config_update)
}
\arguments{
  \item{config_filename}{Path to a pipeline configuration
  file}

  \item{config_update}{A list of named value pairs that
  will update the config parameters}
}
\value{
  The path to the NGS output directory.
}
\description{
  Run the NGS analysis pipeline from a configuration file
}
\details{
  This is the launcher function for all pipeline runs. It
  will do some preprocessing steps, then aligns the reads,
  counts overlap with genomic Features such as genes, exons
  etc and applies a variant caller.
}
\author{
  Jens Reeder, Gregoire Pau
}
\seealso{
  runPipeline
}
\keyword{internal}

