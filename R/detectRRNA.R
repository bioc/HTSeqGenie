##' Returns a named vector indicating if a read ID has rRNA contamination or not 
##'
##' Given a genome and fastq data, each read in the fastq data is
##' aligned against the rRNA sequences for that genome
##' @title Detect rRNA Contamination in Reads
##' @param lreads A list of ShortReadQ objects
##' @param remove_tmp_dir boolean indicating whether or not to delete temp directory of gsnap results
##' @param save_dir Save directory
##' @return a named logical vector indicating if a read has rRNA contamination
##' @author Cory Barr
##' @export
detectRRNA <- function (lreads, remove_tmp_dir=TRUE, save_dir=NULL) {
  ## init
  if (missing(save_dir)) save_dir <- file.path(getConfig("save_dir"))
  setUpDirs(save_dir, overwrite="overwrite")
  loginfo("detectRRNA.R/detectRRNA: searching for rRNA contamination...")
 
  ## force lreads to be a list
  if (!is.list(lreads)) lreads = list(lreads)
  
  ## get config parameters
  paired_ends <- getConfig.logical("paired_ends")
  rRNADb <- getConfig("detectRRNA.rrna_genome")
  
  ## misc
  tmp_rRNA_contam_dir <- file.path(save_dir, "bams", "rRNA_contam")
  tmp_rRNA_contam_dir <- makeDir(tmp_rRNA_contam_dir, overwrite="erase")
  
  ## write forward reads for alignement
  filepaths <- writeFastQFiles(lreads, dir=tmp_rRNA_contam_dir, filename1="input1.fastq", filename2="input2.fastq")

  if(length(filepaths) == 1){
    file2 <- NULL
  } else{
    file2 <- filepaths[[2]]
  }
  
  ## now run gsnap and grep out all ids of rRNA reads
  contam_ids = getRRNAIds(
    file1         = filepaths[[1]],
    file2         = file2,
    tmp_dir       = tmp_rRNA_contam_dir,
    rRNADb        = rRNADb)

  if (remove_tmp_dir) unlink(tmp_rRNA_contam_dir, recursive=TRUE)

  all_ids <- as.character(id(lreads[[1]]))
  
  ## if paired_ends, need to remove the /[12] pair indicator off ids
  if(paired_ends)
    all_ids <- sub("/[12]$", "", all_ids)

  all_ids <- sub(" .*", "", all_ids)
  
  is_contaminated <- all_ids %in% contam_ids
  names(is_contaminated) <- all_ids

  if(sum(is_contaminated)!=length(contam_ids)) {
    stop("ERROR: rRNA contam removal failed.")
  }

  loginfo(paste("detectRRNA.R/detectRRNA: contaminated fraction=", mean(is_contaminated)))
  loginfo("detectRRNA.R/detectRRNA: done")
  return(is_contaminated)
}

##' Detect reads that look like rRNA
##'
##' @param file1 FastQ file of forward reads
##' @param file2 FastQ of reverse reads in paired-end sequencing, NULL otherwise
##' @param tmp_dir temporary directory used for storing the gsnap results
##' @param rRNADb Name of the rRNA sequence database. Must exist in the gsnap genome directory
##' @return IDs of reads flagged as rRNA
##' @export
getRRNAIds <- function (file1,
                        file2=NULL,
                        tmp_dir,
                        rRNADb) {
  ## get config parameters
  path.gsnap_genomes <- getConfig("path.gsnap_genomes")
  max_mismatches <- getConfig.integer("alignReads.max_mismatches")
  prepend_str <- getConfig('prepend_str')

  ## prepare gsnapParams
  gsnapParams <- paste("-t", 0, ## zero threads means main and IO runs in same thread
                       "-D", path.gsnap_genomes,
                       "-n 1",
                       "-d", rRNADb,
                       "-A sam",
                       "-B 2", ## 2 = memory mapped (to share the memory map among parallel gsnap processes)
                       "--pairmax-rna=200000")
  if (!is.null(max_mismatches)) gsnapParams <- paste(gsnapParams, "-m", max_mismatches)

  ## align reads
  outbams <- wrapGsnap(file1, file2, gsnapParams, tmp_dir, prepend_str, remove.nomapping=TRUE)
  
  ## get IDs of rRNA contaminated reads
  sb_param <- ScanBamParam(what=c("qname"))
  contam_ids <- lapply(outbams,
                       function(x) {
                         scanBam(x,
                                 param=sb_param,
                                 tag="MD")})
  contam_ids <- unlist(contam_ids, use.names=FALSE)
  ## dedup because of multimappers
  contam_ids <- contam_ids[!duplicated(contam_ids)]

  return(contam_ids)
}
