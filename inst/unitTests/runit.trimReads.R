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

  ## check that the trimmed reads are the same ones produced by substr(...)
  a <- sread(reads)
  b <- sread(treads)
  checkTrue(all(substr(a, 1, nchar(b))==b))
}

test.truncateReads.trim5 <- function() {
  ##  make some random ShortReads
  reads <- HTSeqGenie:::makeRandomSreads(10, 100)
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(10, 80))
  treads <- truncateReads(reads, 50, 10)

  ## check there are 20 reads of 50 nucleotides
  checkEquals(length(treads), 20)
  checkTrue(all(width(treads)==50))

  ## check that the 10 first nucleotides were trimmed
  a <- sread(reads)
  b <- sread(treads)
  c <- paste0(substr(a, 1, 10), b)
  checkTrue(all(c==substr(a, 1, nchar(c))))

  ##  make some random ShortReads, some shorter than the desired length
  reads <- HTSeqGenie:::makeRandomSreads(10, 100)
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(10, 50))
  treads <- truncateReads(reads, 75, 10)

  ## basic checks
  checkEquals(length(treads), 20)
  checkTrue(all(width(treads)<=width(reads)))
  checkTrue(all(width(treads)<=75))
  
  ## check that the 10 first nucleotides were trimmed
  a <- sread(reads)
  b <- sread(treads)
  c <- paste0(substr(a, 1, 10), b)
  checkTrue(all(c==substr(a, 1, nchar(c))))

  ##  make some random ShortReads, some shorter than the desired trim5
  reads <- HTSeqGenie:::makeRandomSreads(10, 100)
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(10, 20))
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(2, 0))
  reads <- append(reads, HTSeqGenie:::makeRandomSreads(2, 1))
  treads <- truncateReads(reads, 60, 30)
  checkTrue(sum(width(treads)==60)==10 && sum(width(treads)==0)==14)
  treads <- truncateReads(reads, 75, 30)
  checkTrue(sum(width(treads)==70)==10 && sum(width(treads)==0)==14)
  treads <- truncateReads(reads, 75, 19)
  checkTrue(sum(width(treads)==75)==10 && sum(width(treads)==1)==10 && sum(width(treads)==0)==4)
  treads <- truncateReads(reads, 1, 0)
  checkTrue(sum(width(treads)==0)==2 && sum(width(treads)==1)==22)

  ##  test with trim5 only
  treads <- truncateReads(reads, trim5=20)
  checkTrue(sum(width(treads)==0)==14 && sum(width(treads)==80)==10)
  treads <- truncateReads(reads, trim5=1)
  checkTrue(all(substring(as.character(sread(reads)), 2)==as.character(sread(treads))))
  checkTrue(sum(width(treads)==0)==4 && sum(width(treads)==19)==10 &&  sum(width(treads)==99)==10)
  treads <- truncateReads(reads)
  checkTrue(all(sread(treads)==sread(reads)))
}
