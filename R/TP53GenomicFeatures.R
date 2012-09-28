##' Build the genomic features of the TP53 demo region
##'
##' Returns a list of genomic features (gene, exons, transcripts) annotating
##' a region of UCSC hg19 sequence centered on the region of the TP53 gene, with 1 Mb flanking
##' sequence on each side. This is intended as a test/demonstration to run the NGS pipeline
##' in conjunction with the \code{LungCancerLines} data package.
##' 
##' @title Demo genomic features around the TP53 gene
##' @return A list of GRanges objects containing the genomic features 
##' @author Gregoire Pau
##' @seealso TP53Genome, buildGenomicFeaturesFromTxDb, runPipeline
##' @export
TP53GenomicFeatures <- function() {
  ## init
  localDir <- "~/.local/share/HTSeqGenie"
  genomicFeatures <- file.path(localDir, "gfeatures-p53Genome.RData")

  ## is genomic_features file already absent?
  if (!file.exists(genomicFeatures)) {
    if (!file.exists(dirname(genomicFeatures))) dir.create(dirname(genomicFeatures), recursive=TRUE)
    gmapR:::checkPackageInstalled("org.Hs.eg.db", required = TRUE)
    gmapR:::checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene", required = TRUE)
    
    ## get p53 region
    gene <- "TP53"
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    roi <- gmapR:::getGeneRoi(txdb, org.Hs.eg.db::org.Hs.eg.db, gene)
      
    ## build genomic features
    genomic_features <- buildGenomicFeaturesFromTxDb(txdb)
    
    ## subset genomic features to the TP53 region
    genomic_features <- lapply(genomic_features, function(gf) {
      suppressWarnings(gmapR:::subsetRegion(gf, roi, gene))
    })
    
    ## save genomic features
    save(genomic_features, file=genomicFeatures)
  }
  genomicFeatures
}
