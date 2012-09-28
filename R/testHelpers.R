##' setup test framework
##' 
##' @param config.filename configuration file
##' @param config.update update list of config values
##' @param testname name of test case
##' @param package name of package
##' @return the created temp directory
##' @export
##' @keywords internal
setupTestFramework <- function(config.filename, config.update=list(), testname="test",
                               package="HTSeqGenie") {

  tp53Genome <- TP53Genome()
  tp53GenomicFeatures <- TP53GenomicFeatures()

  template_config <- file.path(tempdir(), paste(testname,"config.txt", sep="-"))
  buildConfig(template_config,
              ## input
              paired_ends=TRUE,
              quality_encoding="illumina1.8",
            
              ## aligner
              path.gsnap_genomes=path(directory(tp53Genome)),
              alignReads.genome=genome(tp53Genome),
              alignReads.additional_parameters="--indel-penalty=1 --novelsplicing=1 --distant-splice-penalty=1",
              
              ## gene model
              path.genomic_features=dirname(tp53GenomicFeatures),
              countGenomicFeatures.gfeatures=basename(tp53GenomicFeatures)
              )            
  
  ## build config.filename paths
  config.filename <- getPackageFile(config.filename, package)

  config.update <- c(template_config = template_config,
                     config.update)
  ## load original configuration file and update
  loadConfig(config.filename)
  updateConfig(config.update)

  ## This test's tmp dir
  save_dir <- tempfile(pattern=paste(testname, ".", sep=""), tmpdir=tempdir())

  ## prepare update config and update paths
  config.update$save_dir <- save_dir
  filename1 <- getPackageFile(gsub("inst/", "", getConfig("input_file")), package)
  config.update$input_file <- filename1
  if (isConfig("input_file2") && getConfig.logical("paired_ends")) {
    filename2 <- getPackageFile(gsub("inst/", "", getConfig("input_file2")), package)
    config.update$input_file2 <- filename2
  }

  ## init pipeline
  initPipelineFromConfig(config.filename, config.update)
  
  return(save_dir)
}

## generate a couple if random ShortReadQ, intended for testing
makeRandomSreads <- function(num, len) {
  rstring <- apply(matrix(sample(DNA_ALPHABET[1:4], num * len, replace = TRUE), nrow = num),
                   1, paste, collapse = "")
  randomStrings <- DNAStringSet(rstring)

  qual <- paste(rep("A", len), collapse = "")
  qualities <- FastqQuality(rep(qual, num))
  ids <- BStringSet(as.character(seq_len(num)))
  
  randomStrings <- ShortReadQ(sread = randomStrings, quality = qualities,
                              id = ids)
  return(randomStrings)
}


## Convenience function to write a set of seqs to random fastq                    
writeFastqFromRaw <- function(seqs, quals, ids, save_dir=tempdir()){
  if(!file.exists(save_dir))
    stop("save_dir ", save_dir," does not exists!")
  
  fastq <- ShortReadQ(sread   = DNAStringSet(seqs),
                      quality = FastqQuality(quals),
                      id      = BStringSet(ids))
    
  file <- paste(tempfile(tmpdir=save_dir), "fastq", sep=".")
  writeFastq(fastq, file, mode="w")
  return(file)
}
