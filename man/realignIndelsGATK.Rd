\name{realignIndelsGATK}
\alias{realignIndelsGATK}
\title{realign indels via GATK}
\usage{
realignIndelsGATK(bam.file)
}
\arguments{
  \item{bam.file}{Path to bam.file}
}
\value{
Path to realigned bam file
}
\description{
Realing indels using the GATK tools RealignerTargetCreator
and IndelRealigner. Requires a GATK compatible genome with
a name matching the alignment genome to be installed in
'path.gatk_genome'
}
\author{
Jens Reeder
}

