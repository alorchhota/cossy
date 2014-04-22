readExpression <- function(gctfile, chipfile=NA){
  gct <- read.table(gctfile, skip=2, header=TRUE, sep="\t", quote="", comment.char="", as.is=TRUE)
  colnames(gct) <- tolower(colnames(gct))
  
  if(!is.na(chipfile)){
    ## column names: "Probe Set ID"  "Gene Symbol"  "Gene Title"
    annotations <- read.table(chipfile, header=TRUE, sep="\t", colClasses=c('character'), comment.char="", quote="");
    colnames(annotations) <- tolower(colnames(annotations))
    
    probeAnnotations <- merge(x=gct[,"name",drop=F], y=annotations, by.x="name", by.y="probe.set.id", all.x=TRUE, all.y=FALSE)
    rownames(probeAnnotations) <- probeAnnotations$name
    kidInputs <- probeAnnotations[gct$name,"gene.symbol"]
  }
  else {
    kidInputs <- gct$description
  }
  
  # handle white spaces, ///, --- in kid accordingly
  kidInputs <- gsub("/{2,}", ";", kidInputs)
  kidInputs <- gsub("\\s", "", kidInputs)
  kidInputs <- gsub("-{2,}", "-", kidInputs)
  kidInputs[is.na(kidInputs)] <- "-"
  kidInputs[kidInputs==""] <- "-"
  kid <- data.frame(kid=kidInputs, stringsAsFactors=F)
  
  return(data.frame(gct, kid, row.names=gct[,"name"]))
}

readClass <- function(clsfile){
  cls <- read.table(clsfile, skip=1, header=F, sep="", quote="",  );
  data.frame(class=cls[!is.na(cls)])  # remove NA to avoid ending spaces
}

preprocessTrainingExpression <- function(expression, qnorm=T, ztrans=T){
  
  ## This file customize the function normalizeQuantiles function from limma package.
  ## customQuantileNormalize function normalizes the dataset keeping one test sample out (for loocv),
  ## saves the normalization parameters,
  ## and then normalize the test sample using the saved parameters.
  
  customQuantileNormalize <- function (A, ties = TRUE) 
  {
    n <- dim(A)
    if (is.null(n)) 
      return(A)
    if (n[2] == 1) 
      return(A)
    O <- S <- array(, n)
    nobs <- rep(n[1], n[2])
    i <- (0:(n[1] - 1))/(n[1] - 1)
    for (j in 1:n[2]) {
      Si <- sort(A[, j], method = "quick", index.return = TRUE)
      nobsj <- length(Si$x)
      if (nobsj < n[1]) {
        nobs[j] <- nobsj
        isna <- is.na(A[, j])
        S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x, 
                         i, ties = "ordered")$y
        O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
      }
      else {
        S[, j] <- Si$x
        O[, j] <- Si$ix
      }
    }
    m <- rowMeans(S)
    for (j in 1:n[2]) {
      if (ties) 
        r <- rank(A[, j])
      if (nobs[j] < n[1]) {
        isna <- is.na(A[, j])
        if (ties) 
          A[!isna, j] <- approx(i, m, (r[!isna] - 1)/(nobs[j] - 
                                                        1), ties = "ordered")$y
        else A[O[!isna, j], j] <- approx(i, m, (0:(nobs[j] - 
                                                     1))/(nobs[j] - 1), ties = "ordered")$y
      }
      else {
        if (ties) 
          A[, j] <- approx(i, m, (r - 1)/(n[1] - 1), ties = "ordered")$y
        else A[O[, j], j] <- m
      }
    }
    
    #A
    return(list(norm=A, quantiles=m))
  }
  
  
  nc <- ncol(expression)
  gdata <- expression[,-c(1,2,nc), drop=F]
  
  quantiles <- NA
  centers <- NA
  scales <- NA
  
  if(qnorm){
    qn <- customQuantileNormalize(gdata)
    gdata <- qn$norm
    #testData[,1] <- customSampleQuantileNormalize(as.vector(testData[,1]), qn$quantiles)
    quantiles <- qn$quantiles
  }
  
  if(ztrans){
    gdata <- t(gdata)
    gdata <- scale(gdata)
    centers <- attr(gdata,"scaled:center")
    scales <- attr(gdata,"scaled:scale")
    gdata <- t(gdata)
    #testData[,1] <- (as.vector(testData[,1]) - centers)/scales
  }
  
  processedExpression <- data.frame(expression[,c(1,2)], gdata[,], expression$kid, stringsAsFactors=F)
  colnames(processedExpression) <- colnames(expression)
  
  return(list(expression=processedExpression, qnorm.quantiles=quantiles, ztrans.centers=centers, ztrans.scales=scales))
}

preprocessTestExpression <- function(preprocessObj, expression){
  customSampleQuantileNormalize <- function(sample, quantiles){
    r <- rank(sample)
    n1 <- length(sample)  # #row in A
    isna <- is.na(sample)
    nobs <- sum(!isna)
    norm <- sample
    i <- (0:(n1 - 1))/(n1 - 1)
    if (nobs < n1){
      norm[!isna,] <- approx(i, quantiles, (r[!isna] - 1)/(nobs - 1), ties = "ordered")$y
    }
    else{
      norm <- approx(i, quantiles, (r - 1)/(n1 - 1), ties = "ordered")$y
    }
    
    return(norm)
  }
  
  nc <- ncol(expression)
  testData <- expression[,-c(1,2,nc), drop=F]
  if(!(length(preprocessObj$qnorm.quantiles)==1 && is.na(preprocessObj$qnorm.quantiles))){
    for(j in 1:ncol(testData)){
      testData[,j]  <- customSampleQuantileNormalize(as.vector(testData[,j]), preprocessObj$qnorm.quantiles)
    }
  }
  
  if( !(length(preprocessObj$ztrans.centers)==1 && is.na(preprocessObj$ztrans.centers) &&  length(preprocessObj$ztrans.scales)==1 && is.na(preprocessObj$ztrans.scales)) ){
    for(j in 1:ncol(testData)){
      testData[,j]  <- (as.vector(testData[,j]) - preprocessObj$ztrans.centers)/preprocessObj$ztrans.scales
    }
  }
  
  processedExpression <- data.frame(expression[,c(1,2)], testData[,], expression$kid, stringsAsFactors=F)
  colnames(processedExpression) <- colnames(expression)
  
  return(list(expression=processedExpression))
  
}

cossy <- function(expression, cls, misset, nmis=5){
  
  raw_expression <- expression
  classes <- cls
  initial_gene_sets <- misset
  k <- nmis
  
  mergeThreshold <- 0.6
  n.threshold.local <- 5
  
  ################ filter to remove expression data with no kegg id (kegg id="-") ################
  hasKid <- !(regexpr("-",raw_expression$kid)>0)
  
  expression.ngenes <- sum(hasKid)
  expression.nsamples <- ncol(raw_expression) - 3   ## Non-sample columns: Name, Description, Kid
  
  if(expression.nsamples < 2){
    stop("Not enough samples.")
  }
  
  if(expression.nsamples != nrow(classes)){
    stop("A mismatch found between the number of samples and corresponding class labels.")
  }
  
  if(length(unique(classes[,1])) != 2){
    stop("There must be exactly two unique class labels.")
  }
  
  ####### set some initial values, that could not be saved before reading the inputs.############
  NoOfClusters <- round(sqrt(expression.nsamples/2))
  uniqueClasses <- sort(as.character(unique(classes[,1])))
  negativeClass <- uniqueClasses[1]
  positiveClass <- uniqueClasses[2]
  
  
  createGeneIdToProbeIdMap <- function(expression){
    ### for all probe id, get the kids (gids) and update a reverse map
    splittedKids <- strsplit(expression$kid,";")
    gid2probeid <- list()
    lapply(1:nrow(expression), function(i){
      probeid <- as.character(expression$name[i])
      storeKidToProbeid <- function(kid){
        if(kid != "-")
          gid2probeid[kid][[1]] <<- c(gid2probeid[kid][[1]], probeid)
        return(NA)
      }
      #debug(storeKidToProbeid)
      lapply(unique(splittedKids[[i]]), storeKidToProbeid)
      
      return(NA)
    })
    
    return(gid2probeid)
  }
  
  ################ process initial_gene_sets. if gene_set does not contain enough probes (<5 probes), discard it. ################################
  ################ this is used only to avoid computation cost due to small size geneset, data is not used in training or testing.  ##############
  #print("processing data... ")
  passedGis <- c()
  expression <- raw_expression[hasKid,]
  #expression <- preprocessDataForLOOCV(expression=raw_expression, fold=-1, normalize=normalize, scale=scale_type, fil=hasKid)
  g2pmap <- createGeneIdToProbeIdMap(expression)
  passedGis <- sapply(initial_gene_sets, function(gis){
    return(length(unique(unlist(g2pmap[gis$genes])))>= n.threshold.local)
  })
  #print(Sys.time())
  initial_gene_sets <- initial_gene_sets[passedGis]
  
  
  trainingClasses <- classes
  trainingExpression <- expression
  
  ################# Storage for top k sets in one fold #########################
  clusterings <- list()
  ents <- c()
  kgis <- list()
  kfreq <- list()
  kLocalFilters <- list()
  
  ###### merge functions #######
  mergeGisSets <- function(allGis, expression){
    
    ### find pairwise intersection score of all gis
    #passedGis <- c()
    
    #cat(paste(Sys.time(), "Merge start: creating initial filters ...\n"))
    
    
    allFilters <- sapply(allGis, function(gis){
      setExpression <- getExpressionData(trainingExpression, gis, onlyValues=FALSE)
      #passedGis <<- c(passedGis,(nrow(setExpression) >= n.threshold.local))
      lfil <- localNGeneFilter(setExpression, classLables=trainingClasses, positiveClass=positiveClass, negativeClass=negativeClass, threshold=n.threshold.local)
      return(lfil)
    })
    
    
    
    ngis <- length(allGis)
    #   if(ngis == 0){
    #     return(list(gis=c(), filters=c()))
    #   }
    if(ngis <= 1){
      return(list(gis=allGis, filters=allFilters))
    }
    
    
    overlap <- matrix(0, nrow=2*ngis, ncol=2*ngis) #max size
    
    #cat(paste(Sys.time(), " Calculating initial intersection scores ...\n"))
    
    t <- 1:ngis
    pairs <- combn(t,2)
    
    intersectionScoreOfPair <- function(p){
      # get the indexes of the current pair
      i1 <- p[1]
      i2 <- p[2]
      
      # get the objects for the current pair
      if(is.na(allGis[[i1]]) || is.na(allGis[[i2]])){
        overlap[i1,i2] <<- 0
      }
      else{
        probes1 <- allFilters[,i1]$probes
        probes2 <- allFilters[,i2]$probes
        
        overlap[i1,i2] <<- intersectionScore(probes1, probes2)
        #overlap[i1,i2] <<- sum(probes1 %in% probes2)/5
      }
      
      return(NA)
    }
    
    apply(pairs, 2, intersectionScoreOfPair)
    
    
    # take the top-scoring pair (>threshold).
    # if more than 1 pair, has highest score, sort by sum of scores
    # merge top pair and add to allGIS.
    # remove individual genes from allGIS. put NULL or NA
    # update the overlap matrix. each row of individual genes = 0, new row and col for new GIS
    # continue this until no genes are merged
    
    deScoreOfPair <- function(p){
      de1 <- sum(allFilters[,p[1]]$scores)
      de2 <- sum(allFilters[,p[2]]$scores)
      score <- de1 + de2
      return(score)
    }
    
    tempCount <- 0
    totalTime <- 0
    
    #cat(paste(Sys.time(), " Mergins genes  ...\n"))
    
    while(TRUE){
      
      #t1 <- Sys.time()
      
      tempCount <- tempCount + 1
      
      maxScore <- max(overlap)
      if(maxScore < mergeThreshold){
        break
      }
      
      #if(tempCount%%50==0){
      #  print(Sys.time())
      #}
      
      
      maxIndexes <- which(overlap==maxScore, arr.ind=T)
      
      if(nrow(maxIndexes) > 1){
        deScores <- as.vector(apply(maxIndexes, 1, deScoreOfPair))
        chooseIndex <- which(deScores==max(deScores), arr.ind=T)
        mergeIndexes <- maxIndexes[chooseIndex[1],]
      }
      else{
        mergeIndexes <- maxIndexes[1,]
      }
      
      # merge top pair and add to allGIS.
      mergedGis <- mergeGis(allGis[[mergeIndexes[1]]],allGis[[mergeIndexes[2]]])
      allGis[[length(allGis)+1]] <- mergedGis
      
      candidateProbes <- unique(c(allFilters[,mergeIndexes[1]]$probes, allFilters[,mergeIndexes[2]]$probes))
      mergedSetExpression <- getExpressionValues(trainingExpression[candidateProbes,])
      
      mergedFilter <- localNGeneFilter(mergedSetExpression, classLables=trainingClasses, positiveClass=positiveClass, negativeClass=negativeClass, threshold=n.threshold.local)
      newAllFilters <- cbind(allFilters,mergedFilter)
      
      rm(allFilters)
      
      allFilters <- newAllFilters
      
      rm(mergedSetExpression, mergedFilter, candidateProbes)
      
      ##### 25 sec #####
      
      # remove individual genes from allGIS and put their overlap score to 0
      allGis[[mergeIndexes[1]]] <- NA
      allGis[[mergeIndexes[2]]] <- NA
      
      # update the overlap matrix. each row of individual genes = 0
      overlap[mergeIndexes[1],] <- 0
      overlap[,mergeIndexes[1]] <- 0
      overlap[mergeIndexes[2],] <- 0
      overlap[,mergeIndexes[2]] <- 0
      
      
      
      # update the overlap matrix. new row and col for new GIS
      #newOverlap <- expandMatrix(overlap, extraRow=1, extraCol=1)
      
      
      
      #rm(overlap)
      #overlap <- newOverlap
      #newGisNo <- nrow(overlap)
      
      ngis <- ngis + 1
      
      # intersectionScoreOfPair saves the score globally
      tempList <- lapply(1:(ngis-1), function(n)intersectionScoreOfPair(c(n,ngis)))
      rm(tempList)
      
      #t2 <- Sys.time()
      #totalTime <- totalTime + (t2-t1)
      
    }
    
    #cat("loop time:", totalTime, "\n" )
    
    #cat(paste(Sys.time(), " creating giss to return...\n"))
    
    #mergedGis, correspondingFilters, expression
    toTake <- !is.na(allGis)
    allMergedGis <- allGis[toTake]
    allMergedFilters <- allFilters[,toTake]
    
    #cat(paste(Sys.time(), " returning giss...\n"))
    return(list(gis=allMergedGis, filters=allMergedFilters))
    
    
    
  }
  
  intersectionScore <- function(probes1, probes2){
    return(sum(probes1 %in% probes2)/length(probes1))
  }
  
  expandMatrix <- function(m, extraRow=1, extraCol=1, defaultValue=0){
    m1 <- matrix(defaultValue,nrow=nrow(m)+extraRow,ncol=ncol(m)+extraCol)
    m1[1:nrow(m),1:ncol(m)] <- m
    return(m1)
  }
  
  
  mergeGis <- function(gis1, gis2){
    genes <- union(gis1$genes, gis2$genes)
    genes <- sort(genes)
    id <- paste(gis1$id, gis2$id, sep="_")
    name <- paste(gis1$name, gis2$name, sep="_")
    gis <- list(id=id, name=name, genes=genes);
    return(gis)
  }
  
  ######## utility function ##########
  getExpressionValues <- function(expression, colIndexes=NULL){
    if(is.null(colIndexes))
      colIndexes <- !colnames(expression) %in% c("name","description", "kid");
    
    return(as.data.frame(expression[, colIndexes]))
  }
  
  entropyCluster <- function(data, noOfClusters, trainingClasses){
    return(hclustCluster(data, noOfClusters, trainingClasses))
  }
  
  hclustCluster <- function(data, noOfClusters, trainingClasses){
    #### cluster #####
    
    #d <- getDistanceMatrix(data=data, method="euclidean")
    d <- dist(data, method="euclidean")
    
    fit <- hclust(d,method="ward")
    clusters <- cutree(fit, k=noOfClusters)
    
    splitData <- split(as.data.frame(data), clusters)
    
    centers <- sapply(1:noOfClusters, function(i){
      center <- apply(splitData[[as.character(i)]], 2, median)
      return(center)
    })
    
    centers <- t(centers)
    #colnames(centers) <- colnames(data)
    
    
    #### calculate entropy ####
    classFrequencyInClusters <- countSamplesInClsters(clusters, trainingClasses)
    ent <- entropy(classFrequencyInClusters)
    
    return(list(centers, classFrequencyInClusters, ent, clusters))
  }
  
  countSamplesInClsters <- function(clusters, classes){
    
    clustersWithClass <- data.frame(cluster=clusters, class=classes)
    
    countClassSamplesInCluster <- function(cluster, class){
      sum(apply(clustersWithClass, 1, function(row) return(trim(row["cluster"])==cluster && trim(row["class"])==class)))
    }
    
    clusters <- unique(clusters)
    pos <- sapply(clusters, function(cl) countClassSamplesInCluster(cl, positiveClass))
    neg <- sapply(clusters, function(cl) countClassSamplesInCluster(cl, negativeClass))
    
    counts <- data.frame(cluster=clusters, pos=pos, neg=neg, row.names=clusters)
    return(counts)
    
  }
  
  entropy <- function(frequencyInCluster){
    # normalize by dividing by total no of samples of each class
    frequencyInCluster <- normalizeClusterFrequency(frequencyInCluster)
    
    #calculate the entropy of each cluster
    clusterEntropy <- function(npos, nneg){
      p <- npos / (npos+nneg)
      if(p==1 || p==0){
        return(0)
      }
      else{
        return(-p*log2(p) - (1-p)*log2(1-p))
      }
    }
    
    totalEntropy <- 0
    denominator <- 2 # for binary classification, sum of all the fractions is 2
    for(i in 1:nrow(frequencyInCluster)){
      row <- frequencyInCluster[i,]
      npos <- row$pos
      nneg <- row$neg
      totalEntropy <- totalEntropy + (clusterEntropy(npos, nneg) * (npos+nneg) / denominator)
    }
    
    return(totalEntropy)
    
  }
  
  normalizeClusterFrequency <- function(frequencyInCluster){
    totalPos <- sum(frequencyInCluster$pos)
    totalNeg <- sum(frequencyInCluster$neg)
    frequencyInCluster$pos <- frequencyInCluster$pos / totalPos
    frequencyInCluster$neg <- frequencyInCluster$neg / totalNeg
    return(frequencyInCluster)
  }
  
  
  trim <- function(str){
    str <- gsub(pattern="^[ \t]*", replacement="", x=str)
    str <- gsub(pattern="*[ \t]*$", replacement="", x=str)
    str
  }
  
  getExpressionData <- function(expression, geneset, onlyValues=TRUE){
    genes <- geneset$genes
    
    #splits <- strsplit(expression$kid, ";")
    #rowIndexes <- sapply(splits, function(kids){
    #  return(any(kids %in% genes))
    #})
    
    rowIndexes <- unique(unlist(g2pmap[genes]))
    
    if(onlyValues){
      colIndexes <- !colnames(expression) %in% c("name","description", "kid");
    }
    else{
      colIndexes <- 1:ncol(expression);
    }
    
    return(expression[rowIndexes, colIndexes, drop=FALSE])
  }
  
  localNGeneFilter <- function(expression, classLables, positiveClass, negativeClass, threshold=5){
    
    expVal <- getExpressionValues(expression)
    
    iqr.test <- function(x,y){
      medx <- median(x)
      medy <- median(y)
      iqrx <- IQR(x)
      iqry <- IQR(y)
      nx <- length(x)
      ny <- length(y)
      
      #### (Med1-Med2)/(S1^2/n1+S2^2/n)
      sx <- iqrx*34.1/50  # 1 sd around mean (total 2sd) contains 68% data.
      sy <- iqry*34.1/50
      score <- (medx-medy)/sqrt(sx^2/nx+sy^2/ny)
      return(score)
    }
    
    ourIqrTest <- function(row){
      posVal <- as.numeric(row[classLables==positiveClass])
      negVal <- as.numeric(row[classLables==negativeClass])
      iqrVal <- iqr.test(x=negVal, y=posVal)
      return(abs(iqrVal))
    }
    
    scores <- apply(expVal, 1, ourIqrTest)
    scores <- as.vector(scores)
    ord <- order(scores,decreasing=TRUE)
    fil <- rep(FALSE, length(scores))
    fil[ord[1:threshold]] <- TRUE
    
    topScores <- scores[ord[1:threshold]]
    topProbes <- rownames(expVal[fil,])
    retObj <- list(index=fil, scores=topScores, probes=topProbes)
    
    return(retObj)
    
  }
  
  getSingleKidPerProbe <- function(expression, probeIds){
    kids <- expression[probeIds,"kid"]
    kidSplits <- strsplit(kids, split=";")
    singleKids <- sapply(kidSplits, function(x)return(x[1]))
    return(singleKids)
  }
  
  
  ################# prepare filter for feature selection #################
  # initialize gene sets
  merged_gene_sets <- mergeGisSets(initial_gene_sets, trainingExpression) 
  gene_sets <- merged_gene_sets$gis
  set_filters <- merged_gene_sets$filters  
  
  
  if(length(gene_sets) < k){
    stop(paste0("Not enough data to form minimum number of miss (", k ,")."))
  }
  
  ####### prepare an index_vector to get the expression values of gene set ##############
  expressionValColIndexes <- !colnames(trainingExpression) %in% c("name","description", "kid");
  
  setIndex <- 0   
  for(set in gene_sets){
    setIndex <- setIndex + 1
    
    lfil <- set_filters[,setIndex]
    setExpression <- getExpressionValues(trainingExpression[lfil$probes,], expressionValColIndexes)
    
    #### at least 5 genes have to be present, otherwise, heatmaps cannot be generated ####
    #### moreover, same many GISs become same in practice and appear top      ############
    if(nrow(setExpression) < n.threshold.local){
      #warnings()
      next
    }
    
    #### cluster ####
    setExpression <- t(setExpression)
    clusterResult <- entropyCluster(setExpression, NoOfClusters, trainingClasses)
    clusterCenters <- clusterResult[[1]]
    classFrequencyInClusters <- clusterResult[[2]]
    ent <- clusterResult[[3]]
    
    #### sort and take top k ####
    clusterings[[length(clusterings)+1]] <- clusterCenters
    ents <- c(ents, ent)
    kgis[[length(kgis)+1]] <- set
    kfreq[[length(kfreq)+1]] <- classFrequencyInClusters
    kLocalFilters[[length(kLocalFilters)+1]] <- lfil
    
    ord <- order(ents)[-(k+1)] #remove (k+1)th item (if there is more than k items)
    clusterings <- clusterings[ord] 
    ents <- ents[ord]
    kgis <- kgis[ord]
    kfreq <- kfreq[ord]
    kLocalFilters <- kLocalFilters[ord]
    
  }
  
  topMis <- list()
  for(gisIndex in 1:k){
    gis <- kgis[[gisIndex]]
    lfil <- kLocalFilters[[gisIndex]]
    freq <- kfreq[[gisIndex]]
    centers <- clusterings[[gisIndex]]
    colnames(centers) <- lfil$probes
    
    
    gisProfiles <- getExpressionData(expression, gis, onlyValues=FALSE)
    allProbes <- rownames(gisProfiles)
    
    geneNames <- getSingleKidPerProbe(expression, allProbes) ## kids are always gene ids
    representativeGeneNames <- getSingleKidPerProbe(expression, lfil$probes) ## kids are always gene ids
    
    mis <- list(objects=gis$gene, 
                pathway=gis$name, 
                probes=allProbes,
                genes=geneNames,
                #profiles=gisProfiles,
                representative.probes=lfil$probes,
                representative.genes=representativeGeneNames,
                centers=centers,
                freq=freq)
    
    
    topMis[[length(topMis)+1]] <- mis
  }
  
  #predparams <- list(mis=kgis, center=clusterings, freq=kfreq, fil=kLocalFilters, g2pmap=g2pmap, getEData=getExpressionData)
  return(list(topmis=topMis, cls=list(pos=positiveClass, neg=negativeClass)))
  
}

predict <- function(cossyobj, expression){
#   g2pmap
#   gis <- kgis[[gisIndex]]
#   centers <- clusterings[[gisIndex]]
#   freq <- kfreq[[gisIndex]]
#   lfil <- kLocalFilters[[gisIndex]]
  
  topmis <- cossyobj$topmis
  cls <- cossyobj$cls
  #getExpressionData <- params$getEData
  
  vote <- function(probes, testExpression, clusterCenters, freqInCluster){
    
    colnames(testExpression) <- tolower(colnames(testExpression))
    rownames(testExpression) <- testExpression[,"name"]
    colIndexes <- !colnames(testExpression) %in% c("name","description", "kid");
    testSamplesExpression <- testExpression[probes, colIndexes, drop=F]
    
    sampleExpression <- (t(testSamplesExpression[,1]))[1,]
    
    # find the colsest cluster index
    distances <- apply(clusterCenters, 1, function(center) euclideanDitance(center, sampleExpression))
    closestCluster <- which.min(distances)
    
    # count no of pos and neg samples in the colses center
    npos <- freqInCluster[as.character(closestCluster), "pos"]
    nneg <- freqInCluster[as.character(closestCluster), "neg"]
    
    # vote
    voteForPos <-  npos / (npos + nneg)
    return(voteForPos)
  }
  
  getExpressionData <- function(expression, genes, onlyValues=TRUE){
    #genes <- geneset$genes
    
    #splits <- strsplit(expression$kid, ";")
    #rowIndexes <- sapply(splits, function(kids){
    #  return(any(kids %in% genes))
    #})
    
    rowIndexes <- unique(unlist(g2pmap[genes]))
    
    if(onlyValues){
      colIndexes <- !colnames(expression) %in% c("name","description", "kid");
    }
    else{
      colIndexes <- 1:ncol(expression);
    }
    
    return(expression[rowIndexes, colIndexes, drop=FALSE])
  }
  
  euclideanDitance <- function(p1, p2){
    tmp <- data.frame(p1=p1, p2=p2)
    return(sqrt(sum((p1-p2)^2)))
  }
  
  ######### vote by top k GIS #########
  voteForPos <- 0
  k <- length(topmis)
  for(gisIndex in 1:k){
    
#     predparams <- list(mis=gis, center=clusterings, freq=kfreq, fil=kLocalFilters, g2pmap=g2pmap, getEData=getExpressionData)
#     return(list(topmis=topMis, predparams=predparams)
    
    probes <- topmis[[gisIndex]]$representative.probes
    centers <- topmis[[gisIndex]]$centers
    freq <- topmis[[gisIndex]]$freq
    
#     centers <- params$center[[gisIndex]]
#     freq <- params$freq[[gisIndex]]
    
    vpForCurrent <-  vote(probes, expression, centers, freq)
    #print(vpForCurrent)
    
    voteForPos <- voteForPos + vpForCurrent
  }
  
  prediction <- NA
  voteForNeg <- k-voteForPos
  if(voteForPos > voteForNeg){
    prediction <- cls$pos
  }  else if(voteForPos < voteForNeg){
    prediction <- cls$neg
  } else {
    stop("Binary voting needed.")
  }

  return(prediction)
  
}
