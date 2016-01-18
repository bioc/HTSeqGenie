test.readRNASeqEnds <- function(){
  ## According to samtools test file has these alns:
  ## highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0 99 chr22 46085591 34 75M = 46085749 233 GAGAAACAGCACCCAGGACTATCTTCCAAAGAGTTCTGGATATCCTAAAGAAATCTTCTCATGCTGTTGAGCTTGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCRG:Z:SRC111111NM:i:0MD:Z:75SM:i:34
  ## highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0 147 chr22 46085749 34 35M3092N40M = 46085591 -233 TGCATAGAGTGTTCTGTGAACCAGAATTCAATCAGGAACTTGGATACGATTGGTGTTGCTGTTGATTTGATTCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCRG:Z:SRC111111NM:i:0MD:Z:75SM:i:34XS:A:+
  bam <- getPackageFile("test-data/test_highqualSpliceInPE2SharedExonWithPE1.bam")
  observed <- readRNASeqEnds(bam, paired_ends=FALSE)

  checkEquals(length(observed), 2, "readRNASeqEnds() returns GRangesList with right number of reads")
  checkEquals(length(observed[[1]]), 1, "... and the first GRange has one range")
  checkEquals(length(observed[[2]]), 2, "... and it splits gapped alns into two ranges")
  checkEquals(start(observed[[1]]), 46085591, "... and has the right start")
  checkEquals(end(observed[[1]]), 46085591 + 75 - 1, "... and has the right end")
  checkEquals(start(observed[[2]]), c(46085749, 46085749+35+3092), "... and has the right start for the second GRange")
}

test.readRNASeqEnds.dupmark <- function() {
  ## this file has 4 reads (2 paired reads), none of them are dupmarked
  bamfile <- getPackageFile("test-data/ndup-unmarked.bam")
  reads <- readRNASeqEnds(bamfile, paired_ends=FALSE)
  checkEquals(length(reads), 4, "readRNASeqEnds() returns the correct number of reads")
  
  ## this is the very same file, with 2 reads marked as duplicated
  bamfile <- getPackageFile("test-data/ndup-marked.bam")
  reads <- readRNASeqEnds(bamfile, paired_ends=FALSE)
  checkEquals(length(reads), 2, "readRNASeqEnds() returns the correct number of reads")
}
