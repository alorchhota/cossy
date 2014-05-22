## This file contains cross validation of cossy.

## load cossy functions
source('cossy.R')
options(warn=2)

## settings
dataset = "prostate3"                       # name of the dataset
network = "pathwayapi_clustersone"    # name of the network
n.fold = 10                           # #fold in cross validation
n.mis = 15                        # #topmis to vote

print(Sys.time())

set.seed(101)

# Firstly, you will need to prepare 2 or 3 files (*.gct, *.cls, and *.chip) depending on your situation. 
# If you include the gene symbols in the 'description' column of *.gct file, you'll not require the *.chip file.
# Besides, I am planning to incorporate the network mis data in the R package. 
# However, this part is not yet complete. So, I am reading this data from a file (kegg.gmt).

gctfile <- paste0("data/", dataset, ".gct")
clsfile <- paste0("data/", dataset, ".cls")
gmtfile <- paste0("data/", network, ".gmt")

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

## prepare the samples in each fold of cross validation
n.samples <- nrow(cls)
randomizedSamples <- sample(1:n.samples, size=n.samples, replace=F)
n.samples.in.fold <- c(rep(floor(n.samples/n.fold), n.fold-n.samples%%n.fold), rep(floor(n.samples/n.fold)+1, n.samples%%n.fold))
fold.start <- cumsum(c(1,n.samples.in.fold))

## perform cross validation
cvresults <- lapply(1:n.fold, function(fold){
  
  ## show current fold number
  print(paste("fold", fold))
  
  ## separate training and test data.
  testSampleNumber <- randomizedSamples[fold.start[fold]:(fold.start[fold+1]-1)]
  trdata <- exdata[,-(2+testSampleNumber),drop=F]
  trclass <- cls[-testSampleNumber,,drop=F]
  
  ## If your data is not pre-processed (i.e. not normalized or not z-transformed), 
  ## you can do so using preprocessTrainingExpression() function. 
  ## You may omit this step if your data is already processed.
  
  preprocobj <- preprocessTrainingExpression(expression=trdata,qnorm=T,ztrans=T)
  trdata <- preprocobj$expression
  
  
  ## build cossy model
  
  csy <- cossy(expression=trdata, cls=trclass, misset=kegg, nmis=n.mis)
  
  
  ## get the top genes
  
  #topgenes <- sapply(csy$topmis, function(mis){return(mis$representative.genes)})
  #topgenes <- t(topgenes)
  #print(topgenes)
  
  
  ## predict test data using the trained cossy model.
  
  # data with only the test sample
  tsdata <- exdata[,c(1,2,2+testSampleNumber, ncol(exdata)),drop=F]  
  tsclass <- cls[testSampleNumber,,drop=F]
  
  # perform the same preprocessing on test data before prediction
  tsdata <- preprocessTestExpression(preprocessObj=preprocobj, expression=tsdata)$expression
  
  # predict test sample and show output
  prediction <- predict(cossyobj=csy,expression=tsdata)
  
  return(prediction)
})


## save all results sorted by original sample index
predictedClasses <- rep("", n.samples)
predictedVotes <- rep(-1, n.samples)
for(fold in 1:n.fold){
  testSampleNumber <- randomizedSamples[fold.start[fold]:(fold.start[fold+1]-1)]
  predictedClasses[testSampleNumber] <- cvresults[[fold]]$cls
  predictedVotes[testSampleNumber] <- cvresults[[fold]]$vote
}

accuracy <- sum(as.character(cls[,1]) == predictedClasses) / n.samples

print(paste0(dataset, "(n.mis=", n.mis, ") Accuracy: ", format(accuracy*100, digits=2, nsmall=2), " %"))

print(Sys.time())