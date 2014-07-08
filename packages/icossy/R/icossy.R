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

#data("pathwayapi")
#data("pathwayapi_json")

if(getRversion() >= "2.15.1")  utils::globalVariables(c("pathwayapi", "pathwayapi_json"))


getMisGraphFromJsonGmt <- function(jsonGmt, misnumber){
  misnumbers <- sapply(jsonGmt, function(mis) mis$misnumber)
  misIndex <- which(misnumbers == misnumber)
  if(length(misIndex)>0){
    return(jsonGmt[[misIndex]]$graph)
  }
  return(NULL)
}

setExpressionStatusInJsonGraph <- function(jsonGraph, genes, estatus){
  if(length(genes)!=length(estatus) || length(genes)==0)
    stop('gene and estatus must be of same non-zero length.')
  
  nodeLabels <- lapply(jsonGraph$data$nodes, function(node) node$label)
  for(gi in 1:length(genes)){
    nodeIndexes <- which(nodeLabels == genes[gi])
    for(ni in nodeIndexes){
      jsonGraph$data$nodes[[ni]]$expression = estatus[gi]
      jsonGraph$data$nodes[[ni]]$color = ifelse(estatus[gi]=="overexpressed", "red", "green")
    }
  }
  
  return(jsonGraph)
}

lazyLoad <- function(dataVar){
  if(!exists(dataVar)){
    if(dataVar == "pathwayapi")
      data("pathwayapi")
    else if(dataVar == "pathwayapi_json")
      data("pathwayapi_json")
    else
      stop(paste('Invalid data to load :', dataVar))
  }
}

getNetworkData <- function(network="pathwayapi"){
  lazyLoad(network)
  switch(network,
         pathwayapi=pathwayapi)
}

getNetworkJsonData <- function(network="pathwayapi"){
  lazyLoad(paste0(network, "_json"))
  
  switch(network,
         pathwayapi=pathwayapi_json)
}

buildIcossyOutput <- function(cossyobj, cls, jsonGmt){
  misObj <- cossyobj$topmis
  posCls <- cossyobj$cls$pos
  negCls <- cossyobj$cls$neg
  classLabels <- as.character(cls[,1])
  
  
  topMISs <- createEmptyMisList()
  for(mi in 1:length(misObj)){
    mis = misObj[[mi]]
    
    # get mis graph
    misGraph = getMisGraphFromJsonGmt(jsonGmt = jsonGmt, misnumber = mis$subnetid)
    
    # set expression status in mis graph
    misExp <- mis$profiles[mis$representative.probes, ]
    expStatus <- getExpressionStatus(misExp, classLabels, posCls, negCls)
    misGraph = setExpressionStatusInJsonGraph(jsonGraph = misGraph, genes = mis$representative.genes, estatus = expStatus)
    
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

getExpressionStatus <- function(expression, classLables, positiveClass, negativeClass){
  
  colIndexes <- !colnames(expression) %in% c("name","description", "kid");
  expVal <- as.data.frame(expression[, colIndexes])
  
  isOverExpressed <- function(row){
    posVal <- as.numeric(row[classLables==positiveClass])
    negVal <- as.numeric(row[classLables==negativeClass])
    medPos <- median(posVal)
    medNeg <- median(negVal)
    return(medPos>medNeg)
  }
  
  expressionStatus <- apply(expVal, 1, isOverExpressed)
  expressionStatus <- tolower(as.character(unlist(expressionStatus)))
  expressionStatus[expressionStatus=="true"] <- "overexpressed"
  expressionStatus[expressionStatus=="false"] <- "underexpressed"
  
  return(expressionStatus)
  
}


icossy <- function(gctfile, chipfile=NA, clsfile, network, nmis, frank=T, qnorm=F, ztrans=F, sig.test="ttest"){
  
  tryCatch({
  
    exdata <- readExpression(gctfile = gctfile, chipfile = chipfile)
    cls <- readClass(clsfile=clsfile)
    miss <- getNetworkData(network)
    jsonGmt <- getNetworkJsonData(network)
    
    # preprocess
    preprocobj <- preprocessTrainingExpression(expression=exdata, frank = frank, qnorm = qnorm, ztrans = ztrans)
    processedData <- preprocobj$expression
    
    # build model
    csy <- cossy.v(expression=processedData, cls=cls, misset=miss, nmis=nmis, one.se=F, mis.consistency=T, sig.test = sig.test)
    
    # create icossy output
    jsonCsyNet <- buildIcossyOutput(cossyobj = csy, cls = cls, jsonGmt = jsonGmt)
    toReturn <- list(status="OK", network=jsonCsyNet)
    
  }, error = function(e) {
    errmsg <- paste0("CossyException: ",conditionMessage(e)) 
    errReturn <- list(status="ERROR", error=errmsg)
    return(errReturn)
  })
  
}

