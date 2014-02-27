\name{excludeVariantsByRegions}
\alias{excludeVariantsByRegions}
\title{Filter variants by regions}
\usage{
excludeVariantsByRegions(variants, mask)
}
\arguments{
  \item{variants}{Variants as Vranges, GRanges or VCF
  object}

  \item{mask}{region to mask, given as GRanges}
}
\value{
The filtered variants
}
\description{
Filter variants by regions
}
\details{
This function can be used to filter variants in a given
region, e.g. low complexity and repeat regions
}
\author{
Jens Reeder
}

