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

  ## test minLength
  reads <- ShortReadQ(DNAStringSet(c("ATGAATAGTC", "ATGAATAGTCACGACG")),
                      quality=FastqQuality(c("AAAAAAAAAA", "AAAAAAAAAAAAAAAA")))
  z <- HTSeqGenie::isAboveQualityThresh(reads, minquality=30, minfrac=0.7, minlength=12)
  checkEquals(z, c(FALSE, TRUE))
}
