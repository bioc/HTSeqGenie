\name{processVariants}
\alias{processVariants}
\title{Process variant callings}
\usage{
  processVariants(variants, reports.dir)
}
\arguments{
  \item{variants}{GrangesList of variants}

  \item{reports.dir}{Directory to store the resulting
  images and html}
}
\value{
  Name of generated html file
}
\description{
  For a given eset of variants, calculates variant depth,
  confusion matrix, transition ratio and the mutation
  matrix
}

