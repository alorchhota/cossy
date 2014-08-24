#############################################################################################
##### This script is to create the output for the internet version of cossy (icossy) ########
##### It creates an intermediate interface to call cossy functions from web server   ########
##### It also formats the output to make it easy for loading in web                  ########
#############################################################################################

###########################################################################
##### This script assumes the following objects in the environment    #####
##### 1) pathwayapi (read gmt file of pathwayapi)                     #####
##### 2) pathwayapi_json (read json gmt file of pathwayapi)           #####
###########################################################################

 
library(jsonlite)
library(cossy)

## LOAD_NETWORK
if(getRversion() >= "2.15.1")  utils::globalVariables(c("pathwayapi", "pathwayapi_json", 
                                                        "kegg", "kegg_json",
                                                        "string", "string_json"))


getMisGraphFromJsonGmt <- function(jsonGmt, misnumber){
  misnumbers <- sapply(jsonGmt, function(mis) mis$misnumber)
  misIndex <- which(misnumbers == misnumber)
  if(length(misIndex)>0){
    return(jsonGmt[[misIndex]]$graph)
  }
  return(NULL)
}

setExpressionStatusInJsonGraph <- function(jsonGraph, genes, estatus, pos.col, neg.col){
  if(length(genes)!=length(estatus) || length(genes)==0)
    stop('gene and estatus must be of same non-zero length.')
  
  nodeLabels <- lapply(jsonGraph$data$nodes, function(node) node$label)
  for(gi in 1:length(genes)){
    nodeIndexes <- which(nodeLabels == genes[gi])
    for(ni in nodeIndexes){
      jsonGraph$data$nodes[[ni]]$fold = estatus[gi]
      jsonGraph$data$nodes[[ni]]$expression = ifelse(estatus[gi]>0, "overexpressed", "underexpressed")
      jsonGraph$data$nodes[[ni]]$color = ifelse(estatus[gi]>0, pos.col, neg.col)
    }
  }
  
  return(jsonGraph)
}

## LOAD_NETWORK
lazyLoad <- function(dataVar){
  if(!exists(dataVar)){
    if(dataVar == "pathwayapi")
      data("pathwayapi")
    else if(dataVar == "pathwayapi_json")
      data("pathwayapi_json")
    else if(dataVar == "kegg")
      data("kegg")
    else if(dataVar == "kegg_json")
      data("kegg_json")
    else if(dataVar == "string")
      data("string")
    else if(dataVar == "string_json")
      data("string_json")
    else
      stop(paste('Invalid data to load :', dataVar))
  }
}

## LOAD_NETWORK
getNetworkData <- function(network="pathwayapi"){
  lazyLoad(network)
  switch(network,
         pathwayapi=pathwayapi,
         kegg=kegg,
         string=string)
}

## LOAD_NETWORK
getNetworkJsonData <- function(network="pathwayapi"){
  lazyLoad(paste0(network, "_json"))
  
  switch(network,
         pathwayapi=pathwayapi_json,
         kegg=kegg_json,
         string=string_json)
}

buildIcossyOutput <- function(cossyobj, cls, jsonGmt){
  misObj <- cossyobj$topmis
  posCls <- cossyobj$cls$pos
  negCls <- cossyobj$cls$neg
  classLabels <- as.character(cls[,1])
  
  
  topMISs <- createEmptyMisList()
  
  ## Firstly color all the nodes available in the dataset and save the graphs.
  misGraphs <- list()
  for(mi in 1:length(misObj)){
    mis = misObj[[mi]]
    
    # get mis graph
    misGraph = getMisGraphFromJsonGmt(jsonGmt = jsonGmt, misnumber = mis$subnetid)
    
    # set expression status in mis graph
    misExp <- mis$profiles[mis$probes, ]
    expFolds <- getExpressionFold(misExp, classLabels, posCls, negCls)
    misGraph = setExpressionStatusInJsonGraph(jsonGraph = misGraph, genes = mis$genes, estatus = expFolds, pos.col="maroon", neg.col="olive")
    
    # save misGraph
    misGraphs[[mi]] <- misGraph
  }
  
  ## Now, color the representative nodes. 
  ## These nodes should be colored after all normal colorings.
  ## Thus, representative nodes are rightly colored even in overlapped networks.
  
  for(mi in 1:length(misObj)){
    mis = misObj[[mi]]
    
    # get mis graph
    misGraph = misGraphs[[mi]]
    
    # set expression status of representative genes in mis graph
    misExp <- mis$profiles[mis$representative.probes, ]
    expFolds <- getExpressionFold(misExp, classLabels, posCls, negCls)
    misGraph = setExpressionStatusInJsonGraph(jsonGraph = misGraph, genes = mis$representative.genes, estatus = expFolds, pos.col="red", neg.col="green")
    
    # add misGraph
    topMISs = misList.addMis(misList = topMISs, misNumber = mi, mis = misGraph)
  }
  
  jsonNet = getJson(topMISs)
  
  return(jsonNet)
}

createEmptyMisList <- function(){
  return(list())
}

misList.addMis <- function(misList, misNumber, mis){
  misNode <- list(misnumber=misNumber, graph=mis)
  misList[[length(misList)+1]] <- misNode
  return(misList)
}

getJson <- function(obj){
  toJSON(obj, auto_unbox = T)
}

getExpressionFold <- function(expression, classLables, positiveClass, negativeClass){
  
  colIndexes <- !colnames(expression) %in% c("name","description", "kid");
  expVal <- as.data.frame(expression[, colIndexes])
  
  getFold <- function(row){
    posVal <- as.numeric(row[classLables==positiveClass])
    negVal <- as.numeric(row[classLables==negativeClass])
    medPos <- median(posVal)
    medNeg <- median(negVal)
    
    multiplier <- ifelse(medPos >= medNeg, 1, -1)
    nom <- max(abs(medPos), abs(medNeg))
    den <- min(abs(medPos), abs(medNeg))
    delta <- 0.00001
    fold = multiplier * nom / (den + delta)
    
    return(fold)
  }
  
  expressionFolds <- apply(expVal, 1, getFold)

  return(expressionFolds)
}


icossy <- function(gctfile, chipfile=NA, clsfile, network, nmis, frank=T, qnorm=F, ztrans=F, sig.test="ttest", mis.consistency=T){
  
  tryCatch({
  
    exdata <- readExpression(gctfile = gctfile, chipfile = chipfile)
    cls <- readClass(clsfile=clsfile)
    miss <- getNetworkData(network)
    jsonGmt <- getNetworkJsonData(network)
    
    # preprocess
    preprocobj <- preprocessTrainingExpression(expression=exdata, frank = frank, qnorm = qnorm, ztrans = ztrans)
    processedData <- preprocobj$expression
    
    # build model
    if(mis.consistency){
      csy <- cossy.v(expression=processedData, cls=cls, misset=miss, nmis=nmis, one.se=F, mis.consistency=T, sig.test = sig.test)
    }
    else{
      csy <- cossy( expression=processedData, cls=cls, misset=miss, nmis=nmis, sig.test = sig.test)
    }
    
    
    # create icossy output
    jsonCsyNet <- buildIcossyOutput(cossyobj = csy, cls = cls, jsonGmt = jsonGmt)
    toReturn <- list(status="OK", network=jsonCsyNet, classes=c(csy$cls$pos, csy$cls$neg))
    
  }, error = function(e) {
    errmsg <- paste0("CossyException: ",conditionMessage(e)) 
    errReturn <- list(status="ERROR", error=errmsg)
    return(errReturn)
  })
  
}

getAllExpressionFolds <- function(gctfile, chipfile=NA, clsfile, frank=T, qnorm=F, ztrans=F){
  
  tryCatch({
    
    exdata <- readExpression(gctfile = gctfile, chipfile = chipfile)
    cls <- readClass(clsfile=clsfile)
    
    # preprocess
    preprocobj <- preprocessTrainingExpression(expression=exdata, frank = frank, qnorm = qnorm, ztrans = ztrans)
    processedData <- preprocobj$expression
    
    # filter probes without any gene map
    processedData <- processedData[processedData$kid!="-", ]
    
    ## get positive class (*** the logic must be same as it is in cossy.R ***)
    uniqueClasses <- sort(as.character(unique(cls[,1])))
    negativeClass <- uniqueClasses[1]
    positiveClass <- uniqueClasses[2]
    
    # calculate fold values of every probe
    expFolds <- getExpressionFold(expression = processedData, classLables = cls, positiveClass = positiveClass, negativeClass = negativeClass)
    
    # return fold values in a data frame
    toReturn <- data.frame(name=processedData$name, gene=processedData$kid, fold=expFolds)
    return(toReturn)    
    
  }, error = function(e) {
    errmsg <- paste0("CossyException: ",conditionMessage(e)) 
    errReturn <- list(status="ERROR", error=errmsg)
    return(errReturn)
  })
  
}

