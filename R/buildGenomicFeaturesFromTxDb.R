##' Build genomic features from a TxDb object
##'
##' @title Build genomic features from a TxDb object
##' @param txdb A TxDb object.
##' @return A list named list of GRanges objects containing the biological entities to account for.
##' @author Gregoire Pau
##' @examples \dontrun{
##'   library("TxDb.Hsapiens.UCSC.hg19.knownGene")
##'   txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
##'   genomic_features <- buildGenomicFeaturesFromTxDb(txdb)
##' }
##' @importMethodsFrom GenomicFeatures transcriptsBy exonsBy cdsBy transcripts exons
##' @importMethodsFrom IRanges psetdiff
##' @export
buildGenomicFeaturesFromTxDb <- function(txdb) {
  ## genes
  genomic_features <- list()
  genomic_features$gene <- genes(txdb)

  ## exonic regions of genes
  gene_exonic <- exonsBy(txdb, "gene")
  genomic_features$gene_exonic <- gene_exonic
  
  ## record the coding regions by gene
  gene_coding <- cdsBy(txdb, "gene")
  genomic_features$gene_coding <- gene_coding
  
  ## counts per exons
  exon <- exons(txdb, columns=c("exon_name", "tx_name"))
  if (all(is.na(values(exon)$exon_name))) { ## name exons if no canonical names are provided
    names(exon) <- paste(seqnames(exon), ":", start(exon), "-", end(exon), sep="")
  }
  else names(exon)<- values(exon)$exon_name
  values(exon)$exon_name <- NULL
  genomic_features$exon <- exon

  ## computing disjoint exonic regions
  genomic_features$exon_disjoint <- generateSingleGeneDERs(txdb)
  
  ## intergenic counts
  unstranded_genes <- genomic_features$gene ##need to call gaps (on unstranded gene coords)
  strand(unstranded_genes) <- '*'
  intergenic_regions <- gaps(unstranded_genes)
  intergenic_regions <- intergenic_regions[strand(intergenic_regions) == '*']
  names(intergenic_regions) <- paste(seqnames(intergenic_regions),
                                     ":",
                                     start(intergenic_regions),
                                     '..',
                                     end(intergenic_regions),
                                     sep='')
  genomic_features$intergenic <- intergenic_regions

  ## introns
  gene_exonic_reduced <- reduce(gene_exonic) ## union overlapping exons
  gene <- transcriptsBy(txdb, "gene")
  intron <- psetdiff(gene[names(gene_exonic_reduced)],
                     gene_exonic_reduced[names(gene_exonic_reduced)]) ## set difference
  gene_name <- rep(names(intron), elementLengths(intron))
  intron <- unlist(intron) ## flatten all introns
  names(intron) <- paste(gene_name, ".", seqnames(intron), ":", start(ranges(intron)), "-", end(ranges(intron)), sep="") ## name introns
  genomic_features$intron <- intron

  genomic_features
}

##' Generate DEXSeq-ready exons
##' 
##' generateSingleGeneDERs() generates exons by:
##' 1) disjoining the whole exon set
##' 2) keeping only the exons of coding regions
##' 3) keeping only the exons that belong to unique genes
##'
##' @title generateSingleGeneDERs
##' @param txdb A transcript DB object
##' @return single gene DERs
##' @importFrom GenomicFeatures exonsBy cdsBy
##' @importMethodsFrom GenomicRanges seqnames start end
##' @importMethodsFrom BiocGenerics unlist
##' @importFrom IRanges findOverlaps disjoin
generateSingleGeneDERs <- function(txdb) {
   ## build coding exons
   tx <- exonsBy(txdb)
   cds <- range(cdsBy(txdb, "tx"))
   tx <- tx[names(cds)]
   tx <- pintersect(tx, cds)
   gene <- exonsBy(txdb, "gene")

   ## disjoin exons
   der <- disjoin(unlist(tx, use.names = FALSE))
   hits <- findOverlaps(der, gene, ignore.strand = TRUE)
   over.single.gene <- countQueryHits(hits) == 1L
   hits.single.gene <- hits[over.single.gene[queryHits(hits)]]
   der.single.gene <- der[queryHits(hits.single.gene)]

   ## annotate
   values(der.single.gene)$geneid <- names(gene)[subjectHits(hits.single.gene)]
   names(der.single.gene) <- paste(seqnames(der.single.gene), ":", start(der.single.gene), "-", end(der.single.gene), sep="")
   der.single.gene
 }
