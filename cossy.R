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
  
}

predict <- function(cossyobj, expression){
  
}
