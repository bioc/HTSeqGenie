\name{gatk}
\alias{gatk}
\title{gatk}
\usage{
gatk(gatk.jar.path = getOption("gatk.path"), method, args)
}
\arguments{
  \item{gatk.jar.path}{Path to the gatk jar file}

  \item{method}{Name of the gatk method, e.g.
  UnifiedGenotyper}

  \item{args}{additional args passed to gatk}
}
\value{
0 for success, stops otherwise
}
\description{
Run a command from the GATK
}
\details{
Execute the GATK jar file using the method specified as
arg. Stops if the command executed fails.
}
\author{
Jens Reeder
}

