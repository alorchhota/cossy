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

fuzzyRankNormalize <- function(expressionData, frank){
  ## frank parameter can be a logical value or a vector of 2 numbers (theta1 & theta2)
  ## for detail about frank, see the following paper:
  ## Lim,K. and Wong,L. (2014) Finding consistent disease subnetworks using PFSNet. Bioinformatics, 30, 189–96.
  
  doFrankNormalization <- F
  theta1 <- -1
  theta2 <- -1
  
  if(is.logical(frank) && length(frank)==1){
    doFrankNormalization <- frank
    theta1 <- 0.05
    theta2 <- 0.15
  }
  else if(is.numeric(frank) && length(frank)==2){
    doFrankNormalization <- T
    theta1 <- frank[1]
    theta2 <- frank[2]
  }
  else{
    stop('frank must be a logical value or a numeric vector of length 2.')
  }
  
  ## return the same data if fuzzy ranking is not done.
  if(!doFrankNormalization){
    return(expressionData)
  }
  
  ##### fuzzy ranking of one sample #####
  # inputs:
  #   sampleExpression : an array of expression values
  # outputs:
  #   rankedExpression : an array of ranked expression values
  
  oneSampleFRank <- function(sampleExpression, theta1, theta2){
    thresholds <- quantile(x=sampleExpression,  probs=c(theta2, 1-theta1))
    
    if(thresholds[1] >= thresholds[2])
      stop('Higher threshold must be greater than lower threshold. Please check theta1 and theta2 values.')
    
    zeroValues <- (sampleExpression <= thresholds[1])
    oneValues <- (sampleExpression >= thresholds[2])
    midValues <- !(zeroValues | oneValues)
    
    sampleExpression[zeroValues] <- 0
    sampleExpression[oneValues] <- 1
    sampleExpression[midValues] <- (sampleExpression[midValues] - thresholds[1])/(thresholds[2]-thresholds[1])
    return(sampleExpression)
    
  }
  
  ## normalize all the samples
  normalizedData <- apply(expressionData, 2, oneSampleFRank, theta1, theta2)
  return(normalizedData)
  
}


preprocessTrainingExpression <- function(expression, frank=F, qnorm=T, ztrans=T){
  
  ## This file customize the function normalizeQuantiles function from limma package.
  ## customQuantileNormalize function normalizes the dataset keeping one test sample out (for loocv),
  ## saves the normalization parameters,
  ## and then normalize the test sample using the saved parameters.
  
  customQuantileNormalize <- function (A, ties = TRUE) 
  {
    n <- dim(A)
    if (is.null(n)) 
      stop('Error in quantile normalization: Matrix dimension is null.')
    if (n[2] == 1) 
      stop('Error in quantile normalization: At least two samples required.')
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
  
  gdata <- fuzzyRankNormalize(gdata, frank)
  
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
  
  return(list(expression=processedExpression, qnorm.quantiles=quantiles, ztrans.centers=centers, ztrans.scales=scales, frank=frank))
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
  
  testData <-  fuzzyRankNormalize(testData, preprocessObj$frank)
  
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
  
  n.threshold.local <- 5
  
  ## if useEntrypyPvalue is True, then pavalue of entropy is calculated
  ## using permutation test (changing class labels 1000 times randomly)
  ## otherwise, normal entropy value is used.
  useEntrypyPvalue <- F
  
  ## if mergeGeneSets is True, then gene sets (MISs) are merged first.
  mergeGeneSets <- F
  mergeThreshold <- 0.6
  
  
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
  
  ###### process and merge functions #######
  processGisSets <- function(allGis, expression, merge=T){
    
    ### find pairwise intersection score of all gis
    
    #cat(paste(Sys.time(), "Merge start: creating initial filters ...\n"))
    
    
    allFilters <- sapply(allGis, function(gis){
      setExpression <- getExpressionData(trainingExpression, gis, onlyValues=FALSE)
      #passedGis <<- c(passedGis,(nrow(setExpression) >= n.threshold.local))
      lfil <- localNGeneFilter(setExpression, classLables=trainingClasses, positiveClass=positiveClass, negativeClass=negativeClass, threshold=n.threshold.local)
      return(lfil)
    })
    
    ## remove GIS having less than 5 probes
    gisLength <- sapply(allFilters['probes',], length)
    passedGis <- gisLength >= n.threshold.local
    allGis <- allGis[passedGis]
    allFilters <- allFilters[,passedGis]
    
    ngis <- length(allGis)
    #   if(ngis == 0){
    #     return(list(gis=c(), filters=c()))
    #   }
    if(ngis <= 1 || merge==F){
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
        # Here, if any deScore is NA/NAN, then chooseIndex will have length 0.
        # Consequently, an error occur to assign allGis[[mergeIndexes[1]]] <- NA
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
    return(sum(probes1 %in% probes2)/n.threshold.local)
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
    if(ncol(data) < n.threshold.local){
      return(list(NA, NA, .Machine$integer.max, NA))
    }
    
    
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
    # normalize by dividing by total no of samples of each class
    classFrequencyInClusters <- normalizeClusterFrequency(classFrequencyInClusters)
    ent <- entropy(classFrequencyInClusters)
    
    #### calculate p-value of entropy ####
    #### by permutating class labels 1000 times ####
    if(useEntrypyPvalue){
      randomEntropy <- function(clusters, trainingClasses){
        # randomize training labels
        nlabels <- nrow(trainingClasses) 
        randomClasses <- trainingClasses[sample(1:nlabels, size=nlabels, replace=F),,drop=F]
        classFrequencyInClusters <- countSamplesInClsters(clusters, randomClasses)
        classFrequencyInClusters <- normalizeClusterFrequency(classFrequencyInClusters)
        ent <- entropy(classFrequencyInClusters)
        return(ent)
      }
      
      randomEntropies <- sapply(1:1000, function(i){
        randomEntropy(clusters, trainingClasses)
      })
      
      den = density(randomEntropies, from=0, to=1)
      pval = sum(den$y[den$x<=ent])/sum(den$y)
      ent = pval
    }
    
    return(list(centers, classFrequencyInClusters, ent, clusters))
  }
  
  countSamplesInClsters_v0 <- function(clusters, classes){
    
    clustersWithClass <- data.frame(cluster=clusters, class=classes)
    
    countSamplesInClsterscountClassSamplesInCluster <- function(cluster, class){
      sum(apply(clustersWithClass, 1, function(row) return(trim(row["cluster"])==cluster && trim(row["class"])==class)))
    }
    
    clusters <- unique(clusters)
    pos <- sapply(clusters, function(cl) countClassSamplesInCluster(cl, positiveClass))
    neg <- sapply(clusters, function(cl) countClassSamplesInCluster(cl, negativeClass))
    
    counts <- data.frame(cluster=clusters, pos=pos, neg=neg, row.names=clusters)
    return(counts)
    
  }
  
  countSamplesInClsters <- function(clusters, classes){
    clustersWithClass <- data.frame(cluster=clusters, class=classes)
    freq = table(clustersWithClass)
    clusters = as.numeric(row.names(freq))
    counts = data.frame(cluster=clusters, pos=as.numeric(freq[clusters, positiveClass]), neg=as.numeric(freq[clusters, negativeClass]), row.names=clusters)
    return(counts)
  }
  
  entropy <- function(frequencyInCluster){
    
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

    denominator <- 2 # for binary classification, sum of all the fractions is 2
    
#     totalEntropy <- 0
#     for(i in 1:nrow(frequencyInCluster)){
#       row <- frequencyInCluster[i,]
#       npos <- row$pos
#       nneg <- row$neg
#       totalEntropy <- totalEntropy + (clusterEntropy(npos, nneg) * (npos+nneg) / denominator)
#     }


    clusterEntropies <- apply(frequencyInCluster, 1, function(row){
      clusterEntropy(row["pos"], row["neg"]) * (row["pos"]+row["neg"]) / denominator
    })
    
    totalEntropy <- sum(clusterEntropies)
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
      iqrx <- IQR(x, na.rm=T)
      iqry <- IQR(y, na.rm=T)
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
      
      ### use ttest pvalue for ranking genes
      #if(length(unique(row)) == 1)
      #  pval <- 1
      #else
      #  pval <- t.test(negVal, posVal)$p.value
      #return(1-abs(pval))
    }
    
    scores <- apply(expVal, 1, ourIqrTest)
    scores <- as.vector(scores)
    ord <- order(scores,decreasing=TRUE,na.last=NA)
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
  processed_gene_sets <- processGisSets(initial_gene_sets, trainingExpression, mergeGeneSets) 
  gene_sets <- processed_gene_sets$gis
  set_filters <- processed_gene_sets$filters
  
  
  if(length(gene_sets) < k){
    stop(paste0("Not enough data to form minimum number of miss (", k ,")."))
  }
  
  ####### prepare an index_vector to get the expression values of gene set ##############
  expressionValColIndexes <- !colnames(trainingExpression) %in% c("name","description", "kid");
  
  setIndex <- 0   
  #print(length(gene_sets))
  for(set in gene_sets){
    setIndex <- setIndex + 1
    #print(setIndex)
    
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
  
  ## process row and column names of expression (make lower case)
  colnames(expression) <- tolower(colnames(expression))
  rownames(expression) <- expression[,"name"]
  nonExpressionColNames <- c("name","description", "kid")
  
  topmis <- cossyobj$topmis
  cls <- cossyobj$cls
  #getExpressionData <- params$getEData
  
  vote <- function(probes, testExpression, clusterCenters, freqInCluster){
    colIndexes <- !colnames(testExpression) %in% nonExpressionColNames;
    testSamplesExpression <- testExpression[probes, colIndexes, drop=F]
    
    sampleExpression <- (t(testSamplesExpression[,1]))[1,]
    
    ## find cluster distances
    distances <- apply(clusterCenters, 1, function(center) euclideanDitance(center, sampleExpression))
    
    ###### weighted vote ######
    # find the colsest cluster index
    closestCluster <- which.min(distances)
    
    # count no of pos and neg samples in the colses center
    npos <- freqInCluster[as.character(closestCluster), "pos"]
    nneg <- freqInCluster[as.character(closestCluster), "neg"]
    
    # vote
    weightedVoteForPos <-  npos / (npos + nneg)
    
    ##### binary vote #####
    ##### find the first unequal closest cluster and vote the majority class in that cluster
    # order clusters according to distanes
    orderedClusters <- order(distances)
    
    # get frequencies in the clusters
    nposs <- freqInCluster[as.character(orderedClusters), "pos"]
    nnegs <- freqInCluster[as.character(orderedClusters), "neg"]
    
    # find which clusters have unequal number of samples in a cluster
    unequalClusters <- which(nposs!=nnegs)
    if(length(unequalClusters)==0)
      stop('Binary voting is not possible: Every cluster contains the same number of samples.')
    firstUnequalCluster <- unequalClusters[1]
    
    # vote to the majority class
    binaryVoteForPos <- ifelse(nposs[firstUnequalCluster] > nnegs[firstUnequalCluster], 1, 0)
    
    return(list(weighted=weightedVoteForPos, binary=binaryVoteForPos))
  }
  
  euclideanDitance <- function(p1, p2){
    tmp <- data.frame(p1=p1, p2=p2)
    return(sqrt(sum((p1-p2)^2)))
  }
  
  #### predict class of all the samples
  k <- length(topmis)
  sampleNames <- setdiff(colnames(expression), nonExpressionColNames)
  predictions <- sapply(sampleNames, function(sampleName){
    colIndexes <- colnames(expression) %in% c(nonExpressionColNames, sampleName);
    testExpression <- expression[, colIndexes, drop=F]
    
    ######### vote by top k GIS #########
    weightedVoteForPos <- 0
    binaryVoteForPos <- 0
    for(gisIndex in 1:k){
      
      probes <- topmis[[gisIndex]]$representative.probes
      centers <- topmis[[gisIndex]]$centers
      freq <- topmis[[gisIndex]]$freq
  
      # vote
      votesForCurrent <-  vote(probes, testExpression, centers, freq)
      weightedVForCurrent <- votesForCurrent$weighted
      binaryVForCurrent <- votesForCurrent$binary
      
      # sum of votes
      weightedVoteForPos <- weightedVoteForPos + weightedVForCurrent
      binaryVoteForPos <- binaryVoteForPos + binaryVForCurrent
    }
    
    prediction <- NA
    finalVoteForPos <- NA
    
    weightedVoteForNeg <- k - weightedVoteForPos
    if(weightedVoteForPos != weightedVoteForNeg){
      prediction <- ifelse(weightedVoteForPos > weightedVoteForNeg, cls$pos, cls$neg)
      finalVoteForPos <- weightedVoteForPos
    }  else {
      ##### Binary voting needed ####
      #print("binary votes!!!")
      binaryVoteForNeg <- k - binaryVoteForPos 
      prediction <- ifelse(binaryVoteForPos > binaryVoteForNeg, cls$pos, cls$neg)
      finalVoteForPos <- binaryVoteForPos
    }
    
    return(c(cls=prediction, vote=finalVoteForPos/k))
    
  })
  
  predictedClasses = setNames(predictions[1,],sampleNames)
  predictedVotes = setNames(as.numeric(predictions[2,]), sampleNames)
  return(list(cls=predictedClasses,vote=predictedVotes))
  
}

#### cossy.v() is same as the cossy() except that 
#### the final n.mis is selected in validation step.
#### 'v' stands for validation.
cossy.v <- function(expression, cls, misset, nmis=seq(1,15,2)){
  
  ## prepare the samples in each fold of cross validation
  n.fold <- 10
  n.samples <- nrow(cls)
  randomizedSamples <- sample(1:n.samples, size=n.samples, replace=F)
  n.samples.in.fold <- c(rep(floor(n.samples/n.fold), n.fold-n.samples%%n.fold), rep(floor(n.samples/n.fold)+1, n.samples%%n.fold))
  fold.start <- cumsum(c(1,n.samples.in.fold))
  
  ## perform cross validation
  cvresults <- lapply(1:n.fold, function(fold){
    
    ## show current fold number
    print(paste("validation fold", fold))
    
    ## separate training and test data.
    validationSampleNumber <- randomizedSamples[fold.start[fold]:(fold.start[fold+1]-1)]
    trdata <- expression[,-(2+validationSampleNumber),drop=F]
    trclass <- cls[-validationSampleNumber,,drop=F]
    
    
    ## build cossy model
    csy <- cossy(expression=trdata, cls=trclass, misset=kegg, nmis=max(nmis))
    
    
    ## predict validation data using the trained cossy model.
    vdata <- expression[,c(1,2,2+validationSampleNumber, ncol(expression)),drop=F]  
    #vclass <- cls[validationSampleNumber,,drop=F]
    
    # predict test sample and show output
    predictions <- lapply(nmis, function(nm){
      csy1 <- csy
      csy1$topmis <- csy1$topmis[1:nm]
      prediction <- predict(cossyobj=csy1,expression=vdata)
    })
    
    return(predictions)
  })
  
  ## save all results sorted by original sample index
  predictedClasses <- list()
  for(nm in nmis){
    predictedClasses[[as.character(nm)]] <- rep("", n.samples)
  }
  
  for(fold in 1:n.fold){
    validationSampleNumber <- randomizedSamples[fold.start[fold]:(fold.start[fold+1]-1)]
    for(i in 1:length(nmis)){
      nm = nmis[i]
      predictedClasses[[as.character(nm)]][validationSampleNumber] <- cvresults[[fold]][[i]]$cls
    }
  }
  
  # get accuracies
  accuracies <- sapply(nmis, function(nm){
    #sum(as.character(cls[,1]) == predictedClasses) / n.samples
    sum(as.character(cls[,1]) == predictedClasses[[as.character(nm)]]) / n.samples
  })
  
  # find #mis corresponding to max accuracy
  misIndex <- which.max(accuracies)
  final.nmis <- nmis[misIndex]
  
  #print(sort(accuracies, decreasing=T))
  
  ## return the model using final.nmis
  return(cossy(expression, cls, misset, final.nmis))
  
}