test.getChrGRangesFromBAM <- function(){

  bam_file <- getPackageFile("test-data/empty.bam")
  observed <- getChrGRangesFromBAM(bam_file)

  checkEquals(length(observed), 84, "getChrGRangesFromBAM() returns correct number of targets")

  first <- head(observed, 1)
  checkEquals(as.character(seqnames(observed))[1], 'chr1',
              "getChrGRangesFromBAM() returns correct target name for first entry")
  checkEquals(end(observed)[1], 249250621, "getChrGRangesFromBAM() returns correct target length")
  checkEquals(as.character(seqnames(observed))[84], 'chrY', "getChrGRangesFromBAM() returns correct target name for last entry")
  checkEquals(end(observed)[84], 59373566, "getChrGRangesFromBAM() returns correct target length")
}


test.isFirstFragment <- function () {

  ## template having multiple segments
  checkTrue(HTSeqGenie:::isFirstFragment(1), "isFirstFragment(): 1 is first fragment")
  ## first fragment
  checkTrue(HTSeqGenie:::isFirstFragment(64), "isFirstFragment(): 64 is first fragment")
  ## last fragment
  checkTrue(!HTSeqGenie:::isFirstFragment(128), "isFirstFragment(): 128 is first fragment")
  ## secondary aln, first fragment
  checkTrue(HTSeqGenie:::isFirstFragment(256+64), "isFirstFragment(): 320 is first fragment")
  ## secondary aln, last segment
  checkTrue(!HTSeqGenie:::isFirstFragment(256+128), "isFirstFragment(): 374 is last fragment") 
}        

