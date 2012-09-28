test.runPipeline <- function() {
  ## get data
  tp53Genome <- TP53Genome()
  tp53GenomicFeatures <- TP53GenomicFeatures()
  fastq1 <- getPackageFile("extdata/H1993_TP53_subset2500_1.fastq.gz", "HTSeqGenie")
  fastq2 <- getPackageFile("extdata/H1993_TP53_subset2500_2.fastq.gz", "HTSeqGenie")

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
  checkTrue(TRUE, "runPipeline() runs without errors")
}
