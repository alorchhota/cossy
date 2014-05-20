## load cossy functions
source('cossy.R')

print(Sys.time())

set.seed(101)

# Firstly, you will need to prepare 2 or 3 files (*.gct, *.cls, and *.chip) depending on your situation. 
# If you include the gene symbols in the 'description' column of *.gct file, you'll not require the *.chip file.
# Besides, I am planning to incorporate the network mis data in the R package. 
# However, this part is not yet complete. So, I am reading this data from a file (kegg.gmt).

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

## You may like to separate training and test data.

testSampleNumber <- c(1,34)
trdata <- exdata[,-(2+testSampleNumber)]
trclass <- cls[-testSampleNumber,,drop=F]

## If your data is not pre-processed (i.e. not normalized or not z-transformed), 
## you can do so using preprocessTrainingExpression() function. 
## You may omit this step if your data is already processed.

preprocobj <- preprocessTrainingExpression(expression=trdata,qnorm=T,ztrans=T)
trdata <- preprocobj$expression


## build cossy model

csy <- cossy(expression=trdata, cls=trclass, misset=kegg, nmis=15)


## get the top genes

topgenes <- sapply(csy$topmis, function(mis){return(mis$representative.genes)})
topgenes <- t(topgenes)
print(topgenes)


## predict test data using the trained cossy model.

# data with only the test sample
tsdata <- exdata[,c(1,2,2+testSampleNumber, ncol(exdata))]  
tsclass <- cls[testSampleNumber,,drop=F]

# perform the same preprocessing on test data before prediction
tsdata <- preprocessTestExpression(preprocessObj=preprocobj, expression=tsdata)$expression

# predict test sample and show output
prediction <- predict(cossyobj=csy,expression=tsdata)
print(paste0("Actual Class: ", tsclass[,'class'], "; Predicted Class: ", prediction$cls))

print(Sys.time())