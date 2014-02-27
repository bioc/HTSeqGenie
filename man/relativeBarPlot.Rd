\name{relativeBarPlot}
\alias{relativeBarPlot}
\title{Make relative bar plots}
\usage{
relativeBarPlot(data, total, labels, title, filename, ylab = "Percent",
  cex.names = 0.9, ymax = 100)
}
\arguments{
  \item{data}{vector of raw, absolute counts}

  \item{total}{number to normalize by, can be vector of
  same length as data}

  \item{labels}{x-axes labels, category labels for data}

  \item{title}{Title of the plot}

  \item{filename}{plots will be saved under [filename].png
  and [filename].pdf}

  \item{ylab}{label of y axis}

  \item{cex.names}{scaling param of lables, passed to plot}

  \item{ymax}{extent of y-axis}
}
\value{
Nothing, creates two files instead
}
\description{
Make relative bar plots
}
\keyword{internal}

