test.markDuplicates_w_outfile<- function(){
  if (checkPicardJar('MarkDuplicates')){
    tmp <- tempfile()
    bam.file <- getPackageFile("test-data/test_picard_dups.bam")
    markDuplicates(bam.file, tmp)
    
    checkTrue(file.exists(tmp), 'markDuplicates() creates file')
    
    reads <- scanBam(tmp, param=ScanBamParam(what=c("flag")))[[1]]
    checkEqualsNumeric(sum(reads$flag > 1024), 2, "markDuplicates() finds one pair of dups.")
  } else {
    DEACTIVATED("Skipped markDuplicates() test")
  }
}

test.markDuplicates <- function(){
  ## check without passing outfile
  ## copy test file to tmp dir so we don't write next to the original test file
  if (checkPicardJar('MarkDuplicates')){
    bam.file <- getPackageFile("test-data/test_picard_dups.bam")
    tmpdir <- HTSeqGenie:::createTmpDir(prefix="test_markDuplicates")
    inp <- file.path(tmpdir, "pic_md.bam")
    file.copy(bam.file, inp)
    
    outfile <- markDuplicates(inp)
    
    checkTrue(file.exists(outfile), 'markDuplicates() creates file')
    reads <- scanBam(tmp, param=ScanBamParam(what=c("flag")))[[1]]
    checkEqualsNumeric(sum(reads$flag > 1024), 2, "markDuplicates finds one pair of dups.") 
  } else {
    DEACTIVATED("Skipped markDuplicates() test")
  }
}
      
