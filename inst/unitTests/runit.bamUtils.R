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

