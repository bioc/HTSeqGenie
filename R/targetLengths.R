##' Plot target length for paired end 
##'
##' Calculate and plot a histogram of mapped target lengths
##' after trimmming of trim/2 of the data points at the lower and upper
##' end of the distribution
##' 
##' @param bamfile Path to a bam file
##' @param save_dir Path to a piplein results dir
##' @param trim Amount of data to be trimmed at the edges
##' @return Target length table and writes two files in {save_dir}/reports/images/TargetLenghts.[pdf|png]"
##' @author Jens Reeder, Melanie Huntley
##' @export
##' @keywords internal
calculateTargetLengths <- function(bamfile, save_dir, trim=0.4){
  ## get isize
  target_len <- scanBam(bamfile, param=ScanBamParam(what="isize"))[[1]]$isize
  target_len <- abs(as.numeric(target_len))
  target_len_summary <- summary(target_len)
  
  ## keep only isize lower than 600
  trimmed_target_len <- target_len[target_len<600]
  trimmed_target_len_summary <- summary(trimmed_target_len)

  plot_hist <- function() {
    hi <- hist(trimmed_target_len, breaks=100, col="deepskyblue4",
               main="Histogram of target lengths",
               xlab="Target length", border="white", xlim=c(0, 600))
               
    abline(v=trimmed_target_len_summary[3], col="grey", lwd=2)
    text(trimmed_target_len_summary[3], max(hi$counts),
         labels=paste("Median =", trimmed_target_len_summary[3], "bp"), pos=4)
  }

  image.dir <- file.path(save_dir, "reports", "images")    
  png(filename=file.path(image.dir,"TargetLengths.png"),
      width=900*1.2, height=480, units="px")
  dev.control(displaylist =  "enable") # allows dev.copy'ing
  plot_hist()
  dev.copy2pdf( width=15*1.2, height=8,file=file.path(image.dir, "TargetLengths.pdf"))
  dev.off()

  ## write target lengths on disk
  df <- data.frame(names=names(target_len_summary), values=c(target_len_summary))
  saveWithID(df, "summary_target_lengths", id=getConfig('prepend_str'),
             save_dir=file.path(save_dir, "results"),
             format="tab")
  
  ## return targetLengthTable
  targetLengthTable <- data.frame(Trimmed_Data_Set=c(trimmed_target_len_summary),
                                  Full_Data_Set=c(target_len_summary))
  rownames(targetLengthTable) <- names(trimmed_target_len_summary)
  return(targetLengthTable)
}

