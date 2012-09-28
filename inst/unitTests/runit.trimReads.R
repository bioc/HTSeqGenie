test.truncateReads <- function() {
  ## make some random ShortReads
  reads <- HTSeqGenie:::makeRandomSreads(10, 100)

  ## check there are 10 reads of 50 nucleotides
  treads <- truncateReads(reads, 50)
  checkEquals(length(treads), 10)
  checkTrue(all(width(treads)==50))

  ## check we can trim with a number of nucleotides higher than the size of the reads
  treads <- truncateReads(reads, 110)
  checkTrue(all(width(treads)<=110))
  
  ## check quality encoding remains the same after truncateReads
  cqreads <- class(quality(reads))
  cqtreads <- class(quality(treads))
  checkEquals(cqreads, cqtreads)

  ## check we can trim variable size reads
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(10, 80))
  treads <- truncateReads(reads, 50)
  checkEquals(length(treads), 20)
  checkTrue(all(width(treads)==50))
}
