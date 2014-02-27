## tests use TP53 genome as fake rRNA genome
## First exon of TP53
##     TTCCATACCTTCCAATTTTTCAAGGGAAGCTAGAAATCTCGATGTTCATGTGACAGCTCTTTTCTTCTTTATTTTTCAAATATTACCTTATTAAGGTATATTTTGCATATTCATATAATTTGCAGTGTATAATAATTTTTTAAATAAATTTACCAAGTTCTGCAACCATCACCATAAATCGGTTCATCAA
## RC: TTGATGAACCGATTTATGGTGATGGTTGCAGAACTTGGTAAATTTATTTAAAAAATTATTATACACTGCAAATTATATGAATATGCAAAATATACCTTAATAAGGTAATATTTGAAAAATAAAGAAGAAAAGAGCTGTCACATGAACATCGAGATTTCTAGCTTCCCTTGAAAAATTGGAAGGTATGGAA

test.detectRRNA <- function() {
  ## setup test framework

  setupTestFramework(config.filename="test-data/test_config_single_end.txt", testname="test.detectRRNA",
                     config.update=list(
                       detectRRNA.rrna_genome = genome( TP53Genome()),
                       detectRRNA.do=TRUE) )

  shortread1 <- ShortReadQ(          DNAStringSet("TTCCATACCTTCCAATTTTTCAAGGGAAGCTAGAAATCTCGATGTTCATGT"),
                           quality = FastqQuality("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"),
                           id      = BStringSet("tp53_substring_1_50#0/1"))
  
  observed <- detectRRNA(list(shortread1, NULL) )

  checkTrue(observed[1], "detectRRNA() detects rRNA in single end mode")
}

test.detectRRNA.paired_end <- function() {
  ## setup test framework
  setupTestFramework(config.filename="test-data/test_config.txt",
                     testname="test.detectRRNA.paired_end",
                     config.update=list(
                       detectRRNA.do=TRUE,
                       detectRRNA.rrna_genome = genome( TP53Genome() )))

  shortread1 <- ShortReadQ(          DNAStringSet("TTCCATACCTTCCAATTTTTCAAGGGAAGCTAGAAATCTCGATGTTCATGT"),
                           quality = FastqQuality("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"),
                           id      = BStringSet("tp53_substring_1_50#0/1"))

  shortread2 <- ShortReadQ(          DNAStringSet("GGCGCCCGCCTGCTCTCGGTGAGCGCACGTCCCGTGCTCCCCTCTGGCGGG"),
                           quality = FastqQuality("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"),
                           id      = BStringSet("tp53_substring_1_50#0/2"))
  
  observed <- detectRRNA(list(shortread1, shortread2))

  checkTrue(observed[[1]], "detectRRNA() detects rRNA in paired end mode")
}

test.getRRNAIds <- function() {
    
  ## make test files
  tmpdir <- HTSeqGenie:::createTmpDir(prefix="test_get_rRNA_ids")
  file1 <- file.path(tmpdir, "1.fastq")
  file2 <- file.path(tmpdir, "2.fastq")
  
  shortread1 <- ShortReadQ(          DNAStringSet("TTCCATACCTTCCAATTTTTCAAGGGAAGCTAGAAATCTCGATGTTCATGT"),
                           quality = FastqQuality("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"),
                           id      = BStringSet("tp53_substring_1_50#0/1"))

  shortread2 <- ShortReadQ(          DNAStringSet("GGCGCCCGCCTGCTCTCGGTGAGCGCACGTCCCGTGCTCCCCTCTGGCGGG"),
                           quality = FastqQuality("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"),
                           id      = BStringSet("tp53_substring_1_50#0/2"))

  writeFastq(shortread1, file=file1, lane='', compress=FALSE)
  writeFastq(shortread2, file=file2, lane='', compress=FALSE)

  ## run the function to test
  ids <- getRRNAIds(file1 = file1,
                    file2 = file2,
                    tmp_dir = tmpdir,
                    rRNADb = genome( TP53Genome() ))

  checkEquals(ids[1], "tp53_substring_1_50#0", "get_rRNA_ids() does detect rRNA seqs")

  unlink(tmpdir, recursive=TRUE)
}

test.getRRNAIds_random <-function() {
  ## random seqs that should not match rRNA
  tmpdir <- HTSeqGenie:::createTmpDir(prefix="test_get_rRNAIds_random")
  file1 <-file.path(tmpdir, "1.fastq")

  random_sreads <- HTSeqGenie:::makeRandomSreads(100,100)
  writeFastq(random_sreads, file=file1, lane='', compress=FALSE)
  ids <- getRRNAIds(file1 = file1,
                    tmp_dir = tmpdir,
                    rRNADb = genome( TP53Genome() ))
  
  checkEquals(length(ids), 0, "get_rRNA_ids() does not find rRNAs in random strings")
  unlink(tmpdir, recursive=TRUE)
}

