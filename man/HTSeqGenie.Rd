\name{HTSeqGenie}

\alias{HTSeqGenie}

\docType{package}

\title{Package overview}

\description{
  The \code{HTSeqGenie} package is a robust and efficient software
  to analyze high-throughput sequencing experiments in a reproducible
  manner. It supports the RNA-Seq and Exome-Seq protocols and provides:
  quality control reporting (using the \code{ShortRead} package), detection of adapter contamination, read
  alignment versus a reference genome (using the \code{gmapR} package), counting reads in genomic
  regions (using the \code{GenomicRanges} package), and read-depth coverage computation. 
}

\section{Package content}{
  To run the pipeline:
  \itemize{
  \item \code{runPipeline}
  }
  
  To access the pipeline output data:
  \itemize{
  \item \code{getTabDataFromFile}
  }
  
  To build the genomic features object:
  \itemize{
  \item \code{buildGenomicFeaturesFromTxDb}
  \item \code{TP53GenomicFeatures}
  }
}

\examples{
  \dontrun{
    ## build genome and genomic features
    tp53Genome <- TP53Genome()
    tp53GenomicFeatures <- TP53GenomicFeatures()
    
    ## get the FASTQ files
    fastq1 <- system.file("extdata/H1993_TP53_subset2500_1.fastq.gz", package="HTSeqGenie")
    fastq2 <- system.file("extdata/H1993_TP53_subset2500_2.fastq.gz", package="HTSeqGenie")
    
    ## run the pipeline
    save_dir <- runPipeline(
        ## input
        input_file=fastq1,
        input_file2=fastq2,
        paired_ends=TRUE,
        quality_encoding="illumina1.8",
        
        ## output
        save_dir="test",
        prepend_str="test",
        overwrite_save_dir="erase",
        
        ## aligner
        path.gsnap_genomes=path(directory(tp53Genome)),
        alignReads.genome=genome(tp53Genome),
        alignReads.additional_parameters="--indel-penalty=1 --novelsplicing=1 --distant-splice-penalty=1",
        
        ## gene model
        path.genomic_features=dirname(tp53GenomicFeatures),
        countGenomicFeatures.gfeatures=basename(tp53GenomicFeatures)
        )
   }
}

\keyword{package}
