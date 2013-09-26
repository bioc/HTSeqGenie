test.isAboveQualityThresh <- function() {
  ## test simple read
  ## Sanger quality (see QualityScore-class): ASCII code - 33 ("A"=32, "<"=27)
  read <- ShortReadQ(DNAStringSet("ATGAATAGTC"), quality=FastqQuality("AAAAAAAAAA")) ## 100% of 32
  z <- HTSeqGenie::isAboveQualityThresh(read, minquality=30, minfrac=0.7)
  checkEquals(z, TRUE)
  
  ## test simple read
  read <- ShortReadQ(DNAStringSet("ATGAATAGTC"), quality=FastqQuality("<AA<AA<A<A")) ## 60% of 32
  z <- HTSeqGenie::isAboveQualityThresh(read, minquality=30, minfrac=0.7) 
  checkEquals(z, FALSE)
  
  ## test simple read
  read <- ShortReadQ(DNAStringSet("ATGCCTAGTC"), quality=FastqQuality("<ABADACA<A")) ## 80% of >=32
  z <- HTSeqGenie::isAboveQualityThresh(read, minquality=30, minfrac=0.7)
  checkEquals(z, TRUE)
}

test.filterByLength <- function() {
  reads <- ShortReadQ(DNAStringSet(c("ATGAATAGTC", "ATGAATAGTCACGACG")),
                      quality=FastqQuality(c("AAAAAAAAAA", "AAAAAAAAAAAAAAAA")))
  z <- filterByLength(c(reads), minlength=12)
  checkEquals(z, c(FALSE, TRUE))

  reads2 <- append(reads[2], reads[1])
  z <- filterByLength(c(reads, reads2), minlength=12, TRUE)
  checkEquals(z, c(FALSE, FALSE)) 
}
  
test.trimTailsByQuality <- function() {
  
  reads <- c(ShortReadQ(DNAStringSet(c("GATCTGATAAATGCACGCATCCCCCCCCGGGGAAGGGGGTCAGCGCCCCGCGGCACTTATTAGACCCAGCATTAC")),
                      quality=FastqQuality(c("??@D;BBDFDDFDBEHF@FEGEEFHIGH:F#############################################"))))
  checkEquals( width(trimTailsByQuality(reads, minqual='!')[[1]]), width(reads[[1]]))


  checkEquals(width(trimTailsByQuality(reads, minqual='F')[[1]]), 28)
  checkEquals(width(trimTailsByQuality(reads, minqual='H')[[1]]), 26)              
  checkEquals(width(trimTailsByQuality(reads, minqual='I')[[1]]), 0)
}
                      
