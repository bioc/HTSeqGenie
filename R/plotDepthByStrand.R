##' Plot Read Depth of Variants by Strand
##'
##' Plot Read Depth of Variants by Strand
##' @title Plot Read Depth of Variants by Strand
##' @param gr A GRanges with additional elementMetadata describing variants as created by HTSeqGenie
##' @param file Where to write the output
##' @return 0
##' @author Cory Barr
##' @export
##' @keywords internal
plotDepthByStrand <- function(gr, file=NULL) {

  transparent <- function (cols, alpha = 0.4) {
    col_mat <- col2rgb(cols, alpha = TRUE)
    col_mat["alpha", ] <- alpha * 255
    result <- rgb(col_mat["red", ],
                  col_mat["green", ],
                  col_mat["blue",],
                  col_mat["alpha", ], maxColorValue = 255)
    names(result) <- names(cols)
    return(result)
  }

  y <- values(gr)$count.pos + values(gr)$count.neg
  x <- values(gr)$count.pos / (values(gr)$count.pos + values(gr)$count.neg)

  if(is.null(file)) {
    plot(x, y, pch=19, col=transparent("blue"))
  } else {
    png(filename=paste(file, "png", sep="."), width = 900, height = 480, units = "px")
    plot(x, y)
    dev.off()      

    pdf(width=15*1.2, height=8, file=paste(file, "pdf", sep="."))
    plot(x, y)
    dev.off()
  }

  return(0)
}
