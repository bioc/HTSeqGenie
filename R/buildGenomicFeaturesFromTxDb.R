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
##' @export
buildGenomicFeaturesFromTxDb <- function(txdb) {
  ## genes
  genomic_features <- list()
  gene <- transcriptsBy(txdb, "gene")
  genomic_features$gene <- gene
  
  ## exonic regions of genes
  gene_exonic <- exonsBy(txdb, "gene")
  genomic_features$gene_exonic <- gene_exonic
  
  ## record the coding regions by gene
  gene_coding <- cdsBy(txdb, "gene")
  genomic_features$gene_coding <- gene_coding
  
  ## get GRanges of transcripts
  transcript <- transcripts(txdb, columns=c("tx_name", "gene_id"))
  names(transcript) <- values(transcript)$tx_name
  values(transcript)$tx_name <- NULL
  genomic_features$transcript <- transcript
  
  ## counts per exons
  exon <- exons(txdb, columns=c("exon_name", "tx_name"))
  if (all(is.na(values(exon)$exon_name))) { ## name exons if no canonical names are provided
    names(exon) <- paste(seqnames(exon), ":", start(exon), "-", end(exon), sep="")
  }
  else names(exon)<- values(exon)$exon_name
  genomic_features$exon <- exon

  ## intergenic counts
  unstranded_genes <- unlist(genomic_features$gene) ##need to call gaps (on unstranded gene coords)
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
  intron <- psetdiff(gene[names(gene_exonic_reduced)], gene_exonic_reduced[names(gene_exonic_reduced)]) ## set difference
  gene_name <- rep(names(intron), elementLengths(intron))
  intron <- unlist(intron) ## flatten all introns
  names(intron) <- paste(gene_name, ".", seqnames(intron), ":", start(ranges(intron)), "-", end(ranges(intron)), sep="") ## name introns
  genomic_features$intron <- intron

  genomic_features
}
