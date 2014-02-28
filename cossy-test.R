gctfile <- "data/cns.gct"
clsfile <- "data/cns.cls"
gmtfile <- "data/kegg.gmt"

readGmtFile <- function(gmt_file_path){
  gmt_file <- file(gmt_file_path, "r")
  lines <- readLines(gmt_file)
  close(gmt_file)
  splittedLines <- strsplit(lines, split='\t')
  lapply(splittedLines, function(set) list(id=set[1], name=set[2], genes=set[-c(1,2)]));
}

exdata <- readExpression(gctfile=gctfile)
cls <- readClass(clsfile=clsfile)
kegg <- readGmtFile(gmtfile)

testSampleNumber <- 2
trdata <- exdata[,-(2+testSampleNumber)]
trclass <- cls[-testSampleNumber,,drop=F]
preprocobj <- preprocessTrainingExpression(expression=trdata,qnorm=T,ztrans=T)
trdata <- preprocobj$expression
csy <- cossy(expression=trdata, cls=trclass, misset=kegg, nmis=5)
 
topgenes <- sapply(csy$topmis, function(mis){return(mis$representative.genes)})
topgenes <- t(topgenes)

tsdata <- exdata[,c(1,2,2+testSampleNumber, ncol(exdata))]
tsdata <- preprocessTestExpression(preprocessObj=preprocobj, expression=tsdata)$expression
tsclass <- cls[testSampleNumber,,drop=F]
prediction <- predict(cossyobj=csy,expression=tsdata)
print(paste0("Actual: ", tsclass, "; Predicted: ", prediction))
