test.isAdapter3.primeEnd <- function() {

  adapter_seq <- "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
  ## Melanie's handcrafted example
   highqualAdapterContamIn3PrimeEnd <-
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
  
  reads <- c( highqualAdapterContamIn3PrimeEnd )
  observed = isAdapter(reads, 17, adapter_seqs = adapter_seq)

  checkTrue(observed[1],
            "isAdapter() does detect adapter in handcrafted example")
}

test.isAdapter <- function() {
  adapter_seq <- "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"
  l           <- nchar(adapter_seq)
  random_dna  <- paste(sample(DNA_ALPHABET[1:4],
                             100, replace = TRUE),
                      collapse = '')

   reads <- c(random_dna,
             adapter_seq,
             #read with full adapater and random
             paste( c( adapter_seq, random_dna),
                   collapse = '' ),

             #read with first part of adapter and random
             paste(c(substr(adapter_seq, 1, 10), random_dna),
                   collapse = '' ),
             
             #read with last part of adapter and random
             paste(c(substr(adapter_seq, l - 10, l), random_dna),
                   collapse = '' )
             )
  
  # simply take 17 as threshold as this seemed to be one of the largest
  # one we ever got from simulations
  observed = isAdapter(reads, 17, adapter_seqs = adapter_seq)

  checkTrue(!observed[1],
            "isAdapter() does not flag random string as Adapter")
  
  checkTrue(observed[2],
            "isAdapter() does detect adapter")

  checkTrue(observed[3],
            "isAdapter() does detect adapter in longer read")

  # Note that this guy is not flagged as adapter contam, and this is
  # indeed what we want here
  checkTrue(!observed[4],
            "isAdapter() does not flag truncated adapter in longer read")

  checkTrue(observed[5],
            "isAdapter() does detect truncated adapter in longer read")
}

test.getAdapterSeqs <- function() {
  adapter <- HTSeqGenie:::getAdapterSeqs(TRUE, TRUE, 1)
  checkEquals(adapter, "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG",
              "getAdapterSeqs() returns wrong adapter (Test1)" )

  adapter <- HTSeqGenie:::getAdapterSeqs(TRUE, TRUE, 2)
  checkEquals(adapter, "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
              "getAdapterSeqs() returns wrong adapter (Test2)" )

  adapter <- HTSeqGenie:::getAdapterSeqs(FALSE, FALSE, 1)
  checkEquals(adapter, "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
              "getAdapterSeqs() returns wrong adapter (Test3)" )

  adapter <- HTSeqGenie:::getAdapterSeqs(FALSE, TRUE, 1)
  checkEquals(adapter, "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG",
              "getAdapterSeqs() returns wrong adapter (Test4)" )

  checkEquals(HTSeqGenie:::getAdapterSeqs(TRUE, TRUE , 1),
              HTSeqGenie:::getAdapterSeqs(TRUE, FALSE, 1),
              "getAdapterSeqs() should ignore forced_paired_end when paired_ends==TRUE")
}
