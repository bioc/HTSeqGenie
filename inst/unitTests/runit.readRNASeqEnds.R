test.readRNASeqEnds <- function(){

  ## According to samtools test file has these alns:
  ## highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0 99 chr22 46085591 34 75M = 46085749 233 GAGAAACAGCACCCAGGACTATCTTCCAAAGAGTTCTGGATATCCTAAAGAAATCTTCTCATGCTGTTGAGCTTGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCRG:Z:SRC111111NM:i:0MD:Z:75SM:i:34
  ## highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0 147 chr22 46085749 34 35M3092N40M = 46085591 -233 TGCATAGAGTGTTCTGTGAACCAGAATTCAATCAGGAACTTGGATACGATTGGTGTTGCTGTTGATTTGATTCTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCRG:Z:SRC111111NM:i:0MD:Z:75SM:i:34XS:A:+
  bam <- getPackageFile("test-data/test.bam")
  observed <- readRNASeqEnds(bam)

  checkEquals(length(observed), 2, "readRNASeqEnds() returns GRangesList with right number of reads")
  checkEquals(length(observed[[1]]), 1, "... and the first GRange has one range")
  checkEquals(length(observed[[2]]), 2, "... and it splits gapped alns into two ranges")
  
  checkEquals(names(observed),
              c("highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0#1",
                "highqualSpliceInPE2SharedExonWithPE1:1:1:1:12#0#2"),
              "... and it returns the right read names")

  checkEquals(start(observed[[1]]), 46085591, "... and has the right start")
  checkEquals(end(observed[[1]]), 46085591 + 75 - 1, "... and has the right end")

  checkEquals(start(observed[[2]]), c(46085749, 46085749+35+3092), "... and has the right start for the second GRange")

  checkTrue(!values(observed[[1]])$junction, "... and first read is no jucntion")
  checkTrue(all(values(observed[[2]])$junction), "... but the second is")
}
