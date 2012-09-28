test.relativeBarPlot<-function(){

  relativeBarPlot(c(1,2,4),
                  total=4,
                  labels=c("x","y","z"),
                  title="Test",
                  filename=file.path(tempdir(), "test.relativeBarPlot"),
                  )
  
  checkTrue(file.exists(file.path(tempdir(), "test.relativeBarPlot.pdf")),
            "relativeBarPlot() makes pdf bar plot")
  checkTrue(file.exists(file.path(tempdir(), "test.relativeBarPlot.png")),
            "relativeBarPlot() makes png bar plot")
}

test.plotDF<-function(){

  df <- data.frame(num  = c(0,1,2),
                   count= c(60,30,10))
  plotDF(df,
         ylab="Test-Y",
         xlab="Test-X",
         filename=file.path(tempdir(), "test.plotDF") )
  
  checkTrue(file.exists(file.path(tempdir(), "test.plotDF.pdf")),
            "plotDF() makes pdf plot")
  checkTrue(file.exists(file.path(tempdir(), "test.plotDF.png")),
            "plotDF() makes png plot")
}
