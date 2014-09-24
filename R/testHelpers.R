##' setup test framework
##' 
##' @param config.filename configuration file
##' @param config.update update list of config values
##' @param testname name of test case
##' @param package name of package
##' @param use.TP53Genome Boolean indicating the use of the TP53 genome as template config  
##' @return the created temp directory
##' @export
setupTestFramework <- function(config.filename, config.update=list(), testname="test",
                               package="HTSeqGenie", use.TP53Genome=TRUE) {

  ## load original configuration file and update
  config.filename <- getPackageFile(config.filename, package)
  loadConfig(config.filename)

  ## Patch in tp53 template
  if(use.TP53Genome){
    tp53_template <- buildTP53GenomeTemplate()
    config.update <- c(template_config = tp53_template,
                       config.update)
  }
  updateConfig(config.update)
  ## prepare update config and update paths
  ## this lets the tests use package data regardless of it being installed or run out of a svn co
  save_dir <- tempfile(pattern=paste(testname, ".", sep=""), tmpdir=tempdir())
  config.update$save_dir <- save_dir
  filename1 <- getPackageFile(gsub("inst/", "", getConfig("input_file")), package)
  config.update$input_file <- filename1
  if (isConfig("input_file2")) {
    filename2 <- getPackageFile(gsub("inst/", "", getConfig("input_file2")), package)
    config.update$input_file2 <- filename2
  }

  ## init pipeline
  initPipelineFromConfig(config.filename, config.update)
  
  return(save_dir)
}

##' Create a tp53 config template
##'
##' @title buildTP53GenomeTemplate
##' @return Path to tp53 template file
##' @author Jens Reeder
##' @keywords internal
##' @importFrom gmapR TP53Genome
##' @importMethodsFrom gmapR directory
buildTP53GenomeTemplate <- function(){

  tp53Genome <- TP53Genome()
  tp53GenomicFeatures <- TP53GenomicFeatures()
  
  template_config <- file.path(tempdir(), "tp53-config.txt")
  
  ## creates tp53-config.txt in tempdir
  buildConfig(template_config,
              template_config='default-config.txt',
              
              ## aligner
              path.gsnap_genomes=path(directory(tp53Genome)),
              alignReads.genome=genome(tp53Genome),
              
              ## gene model
              path.genomic_features=dirname(tp53GenomicFeatures),
              countGenomicFeatures.gfeatures=basename(tp53GenomicFeatures)
              )
  return(template_config)
}

##' create fasta genome file of TP53 genome
##' 
##' @title buildTP53FastaGenome
##' @return Path to tp53 genome directory
##' @author Jens Reeder
##' @importFrom rtracklayer export
##' @importFrom gmapR TP53Genome
##' @keywords internal
##' @export
buildTP53FastaGenome <- function() {
  tp53 <- TP53Genome()
  genome.name <- genome(tp53)
  
  tp53seq <- DNAStringSet(getSeq(tp53))
  names(tp53seq) = "TP53"
  export(tp53seq, file.path(tempdir(),paste0(genome.name,".fa")), format="fasta")

  return(tempdir())
}

##TODO wait for GeneGenome to pop up in gmapR and then merge this functionality with TP53
buildAnyFastaGenome <- function(genes) {
  gnome <- AnyGenome(genes)
  genome.name <- genome(gnome)
  
  tp53seq <- DNAStringSet(getSeq(gnome))
  names(tp53seq) = genes
  fasta.name <- file.path(tempdir(),paste0(genome.name,".fa"))
  export(tp53seq, fasta.name, format="fasta")
  indexFa(fasta.name)
  return(tempdir())
}

## copied and generalized from gmapR, backported into gmapR as GeneGenome
AnyGenome <- function (genes = c("TP53")){
  
  genomeName <- gmapR:::geneGenomeName(paste(genes, collapse="_"))
 
  if (genomeName %in% genome(GmapGenomeDirectory(create = TRUE))) {
    GmapGenome(genomeName)
  }
  else {
    gmapR:::checkPackageInstalled("TxDb.Hsapiens.UCSC.hg19.knownGene")
    gmapR:::checkPackageInstalled("BSgenome.Hsapiens.UCSC.hg19")
    gmapR:::checkPackageInstalled("org.Hs.eg.db")
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    roi <- gmapR:::getGeneRoi(txdb, org.Hs.eg.db::org.Hs.eg.db, genes)
    strand(roi) <- "+"
    p53Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                     roi, as.character = FALSE)
    names(p53Seq) <- genes
    genome <- GmapGenome(genome = p53Seq, name = genomeName,
                         create = TRUE)
    exons <- gmapR:::subsetRegion(exonsBy(txdb), roi, genes)
    spliceSites(genome, "knownGene") <- exons
    
    genome
  }
}


##' Generate a couple if random ShortReadQ, intended for testing
##'
##' @title Generate a couple if random ShortReadQ, intended for testing
##' @param num an integer
##' @param len an integer
##' @return a DNAStringSet
##' @author Gregoire Pau
##' @keywords internal
##' @importMethodsFrom ShortRead FastqQuality ShortReadQ
makeRandomSreads <- function(num, len) {
  rstring <- apply(matrix(sample(DNA_ALPHABET[1:4], num * len, replace = TRUE), nrow = num),
                   1, paste, collapse = "")
  randomStrings <- DNAStringSet(rstring)

  qual <- paste(rep("A", len), collapse = "")
  qualities <- FastqQuality(rep(qual, num))
  ids <- BStringSet(as.character(seq_len(num)))
  
  randomStrings <- ShortReadQ(sread = randomStrings, quality = qualities,
                              id = ids)
  return(randomStrings)
}

##' Hashing function for vector
##'
##' @title Hashing function for vector
##' @param x A vector
##' @return A numeric
##' @author Gregoire Pau
##' @export
hashVector <- function(x) sum(x*(1:length(x)))

##' Hashing function for coverage
##'
##' @title Hashing function for coverage
##' @param cov A SimpleRleList object
##' @return A numeric
##' @author Gregoire Pau
##' @importMethodsFrom BiocGenerics sapply
##' @export
hashCoverage <- function(cov) hashVector(sapply(cov, function(c) sum(runLength(c)*runValue(c))))

##' Hashing function for variants
##'
##' @title Hashing function for variants
##' @param var A GRanges object
##' @return A numeric
##' @author Gregoire Pau
##' @export
hashVariants <- function(var) hashVector(as.numeric(start(ranges(var))))
