test.buildCountsGRangesList <- function () {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  observed <- buildGenomicFeaturesFromTxDb(txdb)
  checkEquals(names(observed),
              c("gene", "gene_exonic", "gene_coding", "exon", "exon_disjoint", "intergenic", "intron"))  

  ## we mostly test that we get back anything
  checkTrue(length(observed$gene)> 20000L)
  checkTrue(length(observed$gene) >= length(observed$gene_coding), "There are more genes in total than coding genes" )
}

test.generateSingleGeneDERs <- function () {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  observed <- HTSeqGenie:::generateSingleGeneDERs(txdb)

  ## none of these test are affirming that the code is actually right...
  checkTrue(isDisjoint(observed),
            "generateSingleGeneDERs() creates a disjoint set")
  checkTrue(length (setdiff(observed, unlist(range(cdsBy(txdb, "tx"))))) == 0,
            "SingleGeneDER is fully containes within cds ranges")
  checkTrue(length (setdiff(observed, unlist(exonsBy(txdb)))) == 0,
            "SingleGeneDER is fully containes within exon ranges")
  checkTrue(length(observed) > 2* length( cdsBy(txdb, "gene")),
            "We have a lot more DERs than genes")

}
