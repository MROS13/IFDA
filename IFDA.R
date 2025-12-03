library(geigen)
library(nloptr)

###########################################################
#Classification Measures

confusion <- function(trueClass, predictedClass) {
  #Confusion Matrix between "true" and "predicted" labels.
  uniqueClass <- sort(unique(na.omit(trueClass)))
  numClass <- length(uniqueClass)
  confusionMatrix <- matrix(0, nrow = numClass, ncol = numClass)
  for (i in seq(1, length(trueClass))) {
    trueIndex <- which(trueClass[i] == uniqueClass)
    predIndex <- which(predictedClass[i] == uniqueClass)
    confusionMatrix[trueIndex, predIndex] <- confusionMatrix[trueIndex, predIndex] + 1 
  }
  rownames(confusionMatrix) <- colnames(confusionMatrix) <- uniqueClass
  return(confusionMatrix)
}

accuracy <- function(confusionMatrix = NULL, trueClass, predictedClass) {
  #Accuracy based on confusion matrix or computed from "true" and "predicted" labels.
  if (!is.null(confusionMatrix)) {acc <- sum(diag(confusionMatrix))/sum(confusionMatrix)}
  else {acc <- sum(trueClass == predictedClass)/length(trueClass)}
  return(acc)
}

precision <- function(confusionMatrix, class) {
  #Precision for a given class from confusion matrix.
  return(as.numeric(confusionMatrix[class, class]/sum(confusionMatrix[,class])))
}

recall <- function(confusionMatrix, class) {
  #Recall for a given class from confusion matrix.
  return(as.numeric(confusionMatrix[class, class]/sum(confusionMatrix[class,])))
}

fGlobal <- function(confusionMatrix, numClass) {
  #Compute F Global Measure from confusion matrix and number of elements per class.
  n <- sum(numClass)
  fRes <- 0
  for (i in seq(1, ncol(confusionMatrix))){
    prec <- precision(confusionMatrix, i)
    rec <- recall(confusionMatrix, i)
    if (is.nan(prec) || is.nan(rec)) {f1 <- 0}
    else if (rec == 0 || prec == 0) {f1 <- 0}
    else {f1 <- 2 * prec * rec/(prec + rec)}
    fRes <- fRes + f1 * numClass[i]/n
  }
  return(as.numeric(fRes))
}

###########################################################
#Data Treatment

dataSplit <- function(C, R, Class, trainPerc, seed = NULL) {
  #Split Dataset of Centers, Ranges and labels into train and test datasets.
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(C)
  trainIndex <- sample(seq(1, n), floor(trainPerc * n), replace = FALSE)
  trainIndex <- sort(trainIndex)
  trainC <- C[trainIndex,]
  trainR <- R[trainIndex,]
  trainClass <- Class[trainIndex]
  testC <- C[-trainIndex,]
  testR <- R[-trainIndex,]
  testClass <- Class[-trainIndex]
  return(list(trainC = as.matrix(trainC), testC = as.matrix(testC), trainR = as.matrix(trainR), 
              testR = as.matrix(testR), trainClass = trainClass, testClass = testClass))
}

###########################################################
#Numerical Classifier

seqCorrect <- function(a, b) {
  #Coerce seq to provide valid sequence.
  if (b < a) {return(NULL)} else {return(seq(a, b))}
}

productM <- function(x, y, M) {
  #Inner product of x and y induced by matrix M.
  as.numeric(t(x) %*% M %*% y)
}

normaliseM <- function(x, M) {
  #Normalise vector x with norm induced by matrix M: ||x||_M = t(x)Mx. If norm is close to 0 return the 0 vector.
  norM <- sqrt(productM(x, x, M))
  if (round(norM, 8) > 0) {return(x/norM)} else {return(0*x)}
}

colnorM <- function(A, M) {
  #Normalise columns of matrix A with norm induced by matrix M.
  apply(A, 2, function(x) {normaliseM(x, M)})
}

deltaK <- function(k) {
  #Coefficients related to the assumption of the micro-data distribution.
  deltas <- c(0, 1/4, 1/8, 1/12, 1/24, 1/36 - dnorm(3)/(6*(2*pnorm(3) - 1)))
  return(deltas[k])
}

mallowSqrdist <- function(delta, cx, cy, rx, ry) {
  #Squared Mallows' distance.
  as.numeric(t(cx - cy) %*% (cx - cy) + delta * t(rx - ry) %*% (rx - ry))
}

fisheRatio <- function(x, delta, Bc, Br, Wc, Wr) {
  #Fisher Ratio evaluated in x.
  as.numeric((t(x) %*% Bc %*% x + delta * abs(t(x)) %*% Br %*% abs(x))/(t(x) %*% Wc %*% x + delta * abs(t(x)) %*% Wr %*% abs(x)))
}

fisheRatioPrime <- function(x, delta, Bc, Br, Wc, Wr) {
  #Derivative of Fisher's Ratio.
  dSign <- diag(sign(x))
  gamma <- as.numeric(t(x) %*% Bc %*% x + delta * t(abs(x)) %*% Br %*% abs(x))
  beta <- as.numeric(t(x) %*% Wc %*% x + delta * t(abs(x)) %*% Wr %*% abs(x))
  f <- gamma/beta
  fPrime <- 2/beta * ((Bc - f*Wc) %*% x + delta * dSign %*% (Br - f*Wr) %*% abs(x))
  return(fPrime)
}

solveNumerical <- function(delta, Bc, Br, Wc, Wr, M, alphaMatrix = matrix(nrow = 1, ncol = 0)) {
  #Numerical Method to find the solution of IFDA, orthogonal to the previous discriminant vectors in columns of alphaMatrix.
  #Uses library "nloptr".
  discIndex <- seqCorrect(1, ncol(alphaMatrix))
  fConstraint <- function(x) {
    res <- c(productM(x, x, M) - 1)
    for (i in discIndex) {
      res <- rbind(res, productM(x, alphaMatrix[,i], M))
    }
    return(res)
  }
  fConstraintPrime <- function(x) {
    res <- c(2 * M %*% x)
    for (i in discIndex) {
      res <- rbind(res, c(M %*% alphaMatrix[,i]))
    }
    return(res)
  }
  capture.output(nlopOpt <- nloptr::nloptr(x0 = rep(1, ncol(Wc)),
                             eval_f = function(x) {-fisheRatio(x, delta, Bc, Br, Wc, Wr)},
                             eval_grad_f = function(x) {-fisheRatioPrime(x, delta, Bc, Br, Wc, Wr)},
                             eval_g_eq = fConstraint,
                             eval_jac_g_eq = fConstraintPrime,
                             opts = list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-12, "maxeval" = 100000)))
  #Stopping condition
  if (nlopOpt$status == -2) {nlopOpt$solution <- NULL}
  return(nlopOpt$solution)
}

mallowsPredict <- function(alphaMatrix, delta, testC, testR, muClassC, muClassR, uniqueClass) {
  #Classify test observations projected with alphaMatrix and the class barycenter with centre muClassC and range muClassR from smallest Mallows' distance.
  nTest <- nrow(testC)
  numClass <- length(uniqueClass)
  distMatrix <- distPercMatrix <- matrix(0, nrow = nTest, ncol = numClass)
  colnames(distMatrix) <- colnames(distPercMatrix) <- uniqueClass
  predictedClass <- c()
  for (obs in seqCorrect(1, nTest)) {
    bestClass <- -1
    bestVal <- Inf
    for (cl in seqCorrect(1, numClass)) {
      mallowsClass <- 0
      for (disc in seqCorrect(1, ncol(alphaMatrix))) {
        alphaDisc <- t(alphaMatrix[,disc])
        mallowsClass <- mallowsClass + 
          mallowSqrdist(delta, alphaDisc %*% testC[obs,], alphaDisc %*% muClassC[cl,], 
                        abs(alphaDisc) %*% testR[obs,], abs(alphaDisc) %*% muClassR[cl,])
      }
      distMatrix[obs, cl] <- sqrt(mallowsClass)
      if (mallowsClass < bestVal) {
        bestVal <- mallowsClass
        bestClass <- uniqueClass[cl]
      }
    }
    predictedClass <- append(predictedClass, bestClass)
    for (cl in seqCorrect(1, numClass)) {
      distPercMatrix[obs, cl] <- sum(distMatrix[obs, -cl])/sum(distMatrix[obs,])
    }
  }
  return(list(predictedClass = predictedClass, distMatrix = distMatrix, distPercMatrix = distPercMatrix))
}

ifda <- function(delta, 
                 trainC = NULL,
                 trainR = NULL,
                 testC = NULL,
                 testR = NULL,
                 trainClass = NULL,
                 testClass = NULL,
                 restriction = NULL,
                 eps = 10^-8,
                 threshold = -1,
                 s = NULL
) {
  #Interval Fisher's Discriminant Analysis method to classify Symbolic Interval Data
  #Inputs: delta - coefficient associated with underlying distribution of the micro-data. delta in [0,1/4];
  #trainC and trainR - train data;
  #testC and testR - test data;
  #trainClass and testClass - labels associated with train and test data;
  #restriction - restriction of orthogonality between the different discriminant vectors: "uncorrelated" - center uncorrelated or "orthogonal" - usual orthogonality;
  #threshold - maximum difference between performance adding a new discriminant vector. -1 to compute all possible discriminant vectors;
  #s - Optional argument to specify number of discriminant vectors.
  
  #Output: alphaMatrix - list of discriminant vectors;
  #predictedClass - predicted classes;
  #Wc, Wr, Bc and Br - within and between matrices of centres and ranges;
  #muClassC and muClassR - centres and ranges of sample class barycentres;
  #measures - performance measures of using the first r discriminant vectors, successively;
  #M - matrix of orthogonality restriction;
  #delta - delta value used;
  #fisherValues - fisherRatio evaluated in discriminant vectors.
  
  ####Initial verification of inputs
  
  #Verify if delta is a proper value.
  if (as.numeric(delta) != delta) stop("delta should be a real number.")
  if (delta < 0 || delta > 1/4) stop("delta should be between 0 and 0.25.")
  
  #Verify orthogonality restriction is properly defined. "uncorrelated" - center uncorrelated or "usual" - usual orthogonality
  if (is.null(restriction) || all(restriction != c("uncorrelated", "usual"))) stop('The restriction must be "uncorrelated" or "usual".')
  
  #Dimension verification for C and R.
  if (!is.null(trainC)) {
    if(!is.matrix(trainC)) stop("trainC should be a matrix.")
    dimTrainC <- dim(trainC)
  } else stop("trainC was not specified.")
  if (!is.null(testC)) {
    if(!is.matrix(testC)) stop("testC should be a matrix.")
    dimTestC <- dim(testC)
  } else stop("testC was not specified.")
  if (!is.null(trainR)) {
    if(!is.matrix(trainR)) stop("trainR should be a matrix.")
    dimTrainR <- dim(trainR)
  } else stop("trainR was not specified.")
  if (!is.null(testR)) {
    if(!is.matrix(testR)) stop("testR should be a matrix.")
    dimTestR <- dim(testR)
  } else stop("testR was not specified.")
  if (!is.null(trainClass)) {
    if(!is.factor(trainClass) && !is.numeric(trainClass) && !is.character(trainClass)) stop("trainClass should be a numeric, factor or character vector.")
    if(!is.null(dim(trainClass))) stop("trainClass should be a vector.")
    nTrainClass <- length(trainClass)
  } else stop("trainClass was not specified.")
  if (!is.null(testClass)) {
    if(!is.factor(testClass) && !is.numeric(testClass) && !is.character(testClass)) stop("testClass should be a numeric, factor or character vector.")
    if(!is.null(dim(testClass))) stop("testClass should be a vector.")
    nTestClass <- length(testClass)
  } else stop("testClass was not specified.")
  if (!all.equal(dimTrainC, dimTrainR)) stop(
    paste0("trainC (", paste(dimTrainC, collapse = ', '), ") and trainR (", paste(dimTrainR, collapse = ', '), ") have different dimensions.")
  )
  if (!all.equal(dimTestC, dimTestR)) stop(
    paste0("testC (", paste(dimTestC, collapse = ', '), ") and testR (", paste(dimTestR, collapse = ', '), ") have different dimensions.")
  )
  if (dimTrainC[2] != dimTestC[2]) stop(
    paste0("trainC (", paste(dimTrainC[2], collapse = ', '), ") and testC (", paste(dimTestC[2], collapse = ', '), ") have different number of variables.")
  )
  if (nTrainClass != dimTrainC[1]) stop(
    paste0("trainC (", paste(dimTrainC[1]), ") and trainClass (", paste(nTrainClass), ") have different number of observations.")
  )
  if (nTestClass != dimTestC[1]) stop(
    paste0("testC (", paste(dimTestC[1]), ") and testClass (", paste(nTestClass), ") have different number of observations.")
  )
  
  #Define overall sample barycenter.
  muC <- colMeans(trainC)
  muR <- colMeans(trainR)
  
  #Define number of variables p.
  p <- dimTrainC[2]
  
  #Verify number of discriminant vectors.
  if (!is.null(s)) {
    if (s != round(s) || s <= 0) stop("s should be a positive integer.")
  }
  
  #Count number of observations in each class.
  uniqueClass <- sort(unique(na.omit(trainClass)))
  numClass <- length(uniqueClass)
  sizeClass <- table(factor(trainClass, levels = uniqueClass))
  
  #Define sample class barycenter in each row.
  muClassC <- muClassR <- matrix(0, nrow = numClass, ncol = p)
  for (j in seqCorrect(1, numClass)) {
    equalClassIndex <- trainClass == uniqueClass[j]
    if (sizeClass[j] == 1) {
      muClassC[j,] <- trainC[equalClassIndex,]
      muClassR[j,] <- trainR[equalClassIndex,]
    } else {
      muClassC[j,] <- colMeans(trainC[equalClassIndex,])
      muClassR[j,] <- colMeans(trainR[equalClassIndex,])
    }
  }
  
  #Define the between matrices B_C and B_R.
  Bc <- Br <- matrix(0, nrow = p, ncol = p)
  for (j in seqCorrect(1, numClass)) {
    Bc <- Bc + sizeClass[j] * (muClassC[j,] - muC) %*% t(muClassC[j,] - muC)
    Br <- Br + sizeClass[j] * (muClassR[j,] - muR) %*% t(muClassR[j,] - muR)
  }
  
  #Define the within matrices W_C and W_R.
  Wc <- Wr <- matrix(0, nrow = p, ncol = p)
  for (j in seqCorrect(1, numClass)) {
    classWc <- classWr <- matrix(0, nrow = p, ncol = p)
    classIndex <- which(trainClass == uniqueClass[j])
    for (h in seqCorrect(1, sizeClass[j])) {
      classWc <- classWc + (trainC[classIndex[h],] - muClassC[j,]) %*% t(trainC[classIndex[h],] - muClassC[j,])
      classWr <- classWr + (trainR[classIndex[h],] - muClassR[j,]) %*% t(trainR[classIndex[h],] - muClassR[j,])
    }
    Wc <- Wc + classWc
    Wr <- Wr + classWr
  }
  
  #Define M according to the restriction.
  if (restriction == "orthogonal") {
    M <- diag(p)
  } else {
    M <- Wc/(sum(sizeClass) - numClass)
    colnames(M) <- NULL
  }
  
  accGlobal <- f1Global <- 0
  measure <- list()
  measure$acc <- measure$f1 <- c()
  
  #Algorithm to find discriminant vectors. For delta = 0, the conventional FDA method is employed.
  if (delta == 0) {
    if (!is.null(s)) {
      numDiscriminants <- min(s, min(numClass - 1, p))
    } else {
      numDiscriminants <- min(numClass - 1, p)
    }
    geig <- geigen::geigenvectors(Bc, Wc, numDiscriminants, eps)
    if (is.null(geig)) stop("The problem is ill-conditioned. No discriminant vector is available.")
    alphaMatrix <- colnorM(geig$vectors, M)
    for (i in seqCorrect(1, ncol(alphaMatrix))){
      alphaMatrixAux <- matrix(alphaMatrix[,seq(1, i)], nrow = p)
      mallowsPredicted <- mallowsPredict(alphaMatrixAux, delta, testC, testR, muClassC, muClassR, uniqueClass)
      confusionMatrix <- confusion(testClass, mallowsPredicted$predictedClass)
      acc <- accuracy(confusionMatrix)
      f1 <- fGlobal(confusionMatrix, sizeClass)
      measure$acc <- c(measure$acc, acc)
      measure$f1 <- c(measure$f1, f1)
      if (((acc - accGlobal) <= threshold || (f1 - f1Global) <= threshold) && i != 1 ) {
        alphaMatrix <- matrix(alphaMatrix[,seq(1, i-1)], nrow = p)
        break
      }
      accGlobal <- acc
      f1Global <- f1
    }
  } else {
    ##Symbolic Problem
    
    #Specify maximum number of discriminant vectors.
    if (!is.null(s)) {
      numDiscriminants <- min(s, p)
    } else {
      numDiscriminants <- p
    }
    
    #Set alphaMatrix to grow in columns.
    alphaMatrix <- matrix(nrow = p, ncol = 0)
    
    for (i in seqCorrect(1, numDiscriminants)) {
      
      alpha <- solveNumerical(delta, Bc, Br, Wc, Wr, M, alphaMatrix)
     
      if (is.null(alpha)) break
      
      alphaMatrix <- cbind(alphaMatrix, c(alpha))
      mallowsPredicted <- mallowsPredict(alphaMatrix, delta, testC, testR, muClassC, muClassR, uniqueClass)
      confusionMatrix <- confusion(testClass, mallowsPredicted$predictedClass)
      acc <- accuracy(confusionMatrix)
      f1 <- fGlobal(confusionMatrix, sizeClass)
      measure$acc <- c(measure$acc, acc)
      measure$f1 <- c(measure$f1, f1)
      if (((acc - accGlobal) <= threshold || (f1 - f1Global) <= threshold) && i != 1) {
        alphaMatrix <- matrix(alphaMatrix[,seq(1, ncol(alphaMatrix) - 1)], nrow = p)
        break
      }
      accGlobal <- acc
      f1Global <- f1
    }
  }
  mallowsPredicted <- mallowsPredict(alphaMatrix, delta, testC, testR, muClassC, muClassR, uniqueClass)
  
  fisherValues <- c()
  for (i in seqCorrect(1, ncol(alphaMatrix))) {
    fisherValues <- c(fisherValues, fisheRatio(alphaMatrix[,i], delta, Bc, Br, Wc, Wr))
  }
  
  return(list(alphaMatrix = alphaMatrix, predictedClass = mallowsPredicted$predictedClass, 
              Bc = Bc, Br = Br, Wc = Wc, Wr = Wr, muClassC = muClassC, muClassR = muClassR, 
              distMatrix = mallowsPredicted$distMatrix, distPercMatrix = mallowsPredicted$distPercMatrix, 
              measures = measure, M = M, delta = delta, fisherValues = fisherValues))
}

