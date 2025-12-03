#Update: 03/12/2025

##Mosaic Plots and Classmaps for distance-based classification
#Raymaekers J., Rousseeuw P., Hubert M.. "Class maps for visualizing classification results" doi: 10.1080/00401706.2021.1927849
#Functions inspired by library "classmap"

#Use trainlDAC with distance matrix of train data to estimate parameters of transformations
#Use testlDAC with parameters from trainlDAC and distance matrix of test data to compute lDAC and farness of test data
#plotMosaic for mosaic plot of stacked bars
#plotClassmap for plot of lDAC Vs Farness
#plotSilhouette for silhouette plot
#plotFarness for plot of farness per observation index

library(ggplot2)
library(cluster)
library(dplyr)
library(grid)

source("IFDA.R")

trainlDAC <- function(distMatrix, trueClass) {
  #Estimate parameters of lDAC and farness from a distance matrix distMatrix and the true labels trueClass
  #Assumes distance matrix has column names
  ##Returns:
  #lDAC and farness of train dataset
  #gfarness - global farness to find outlier observation to every class
  #farscales and funcsList to use in testlDAC
  
  uniqueClass <- sort(unique(trueClass))
  numClass <- length(uniqueClass)
  lDACValues <- c()
  #Compute lDAC
  for (i in seq(1, nrow(distMatrix))) {
    distTrue <- distMatrix[i, trueClass[i]]
    distPredicted <- min(distMatrix[i, -which(colnames(distMatrix) == trueClass[i])])
    lDACValues <-  c(lDACValues, exp(-distPredicted/2)/(exp(-distPredicted/2) + exp(-distTrue/2)) )
  }
  diss <- distMatrix
  farnDiss <- c()
  for (i in seq(1, nrow(diss))) { #Vector with distances of each obs to correct label
    farnDiss <- c(farnDiss, diss[i, which(trueClass[i] == uniqueClass)])
  }
  omedian <- median(farnDiss) #Median of distances from each object to true class
  farscales <- rep(0, numClass)
  #Standardise distances using medians
  for (g in seq(1, numClass)) {
    yClass <- which(trueClass == uniqueClass[g])
    medClass <- median(diss[yClass, g])
    farscales[g] <- medClass/omedian
    diss[,g] <- diss[,g]/farscales[g]
  }
  #Yeo-Johnson transform
  YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
    lowInd <- which(y < 0)
    highInd <- which(y >= 0)
    if (lambda != 0) {
      y[highInd] <- ((1 + y[highInd])^(lambda) - 1)/lambda
    }
    else {
      y[highInd] <- log(1 + y[highInd])
    }
    if (lambda != 2) {
      y[lowInd] <- -((1 - y[lowInd])^(2 - lambda) - 1)/(2 - lambda)
    }
    else {
      y[lowInd] <- -log(1 - y[lowInd])
    }
    if (stdToo) {
      if (length(y) > 1) {
        locScale <- cellWise::estLocScale(matrix(y, ncol = 1), type = "hubhub")
        zt <- (y - locScale$loc) / locScale$scale
      }
      else {
        zt <- y
      }
    }
    else {
      zt <- NULL
    }
    return(list(yt = y, zt = zt))
  }
  farnessFinal <- rep(0, ncol(diss))
  funcsList <- list()
  probsMatrix <- matrix(0, nrow = nrow(diss), ncol = numClass)
  for (g in seq(1, numClass)) {
    origfarness <- diss[,g]
    classInd <- which(trueClass == uniqueClass[g])
    farness <- diss[classInd, g]
    indnz  <- which(farness > 1e-10)
    farnz  <- farness[indnz]
    farloc <- median(farnz, na.rm = TRUE)
    farsca <- mad(farnz, na.rm = TRUE)
    if (farsca < 1e-10) {farsca <- sd(farnz, na.rm = TRUE)}
    sfar <- scale(farnz, center = farloc, scale = farsca)
    YJout  <- cellWise::transfo(X = sfar, robust = TRUE, standardize = FALSE, checkPars = list(silent = TRUE), quant = 0.90)
    xt      <- YJout$Y
    tfarloc <- median(xt, na.rm = TRUE)
    tfarsca <- mad(xt, na.rm = TRUE)
    lambda  <- YJout$lambdahats
    #Now transform origfarness (to be used in fig):
    origIndnz <- which(origfarness > 1e-10)
    origFarnz <- origfarness[origIndnz]
    zt        <- scale(YJ(scale(origFarnz, farloc, farsca), lambda = lambda, stdToo = FALSE)$yt, tfarloc, tfarsca)
    probs     <- rep(0, length(origfarness))
    probs[origIndnz] <- pnorm(zt)
    farnessFinal[classInd] <- probs[classInd]
    tfunc <- function(farloc, farsca, lambda, tfarloc, tfarsca) {
      force(farloc); force(farsca); force(lambda); force(tfarloc); force(tfarsca)
      tfuncAux <- function(qs) {
        YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
          lowInd <- which(y < 0)
          highInd <- which(y >= 0)
          if (lambda != 0) {
            y[highInd] <- ((1 + y[highInd])^(lambda) - 1) / lambda
          }
          else {
            y[highInd] <- log(1 + y[highInd])
          }
          if (lambda != 2) {
            y[lowInd] <- -((1 - y[lowInd])^(2 - lambda) - 1) / (2 - lambda)
          }
          else {
            y[lowInd] <- -log(1 - y[lowInd])
          }
          if (stdToo) {
            if (length(y) > 1) {
              locScale <- cellWise::estLocScale(matrix(y, ncol = 1), type = "hubhub")
              zt <- (y - locScale$loc) / locScale$scale
            }
            else {
              zt <- y
            }
          }
          else {
            zt <- NULL
          }
          return(list(yt = y, zt = zt))
        }
        qsIndnz <- which(qs > 1e-10)
        qsnz <- qs[qsIndnz]
        zn <- scale(YJ(scale(qsnz, farloc, farsca), lambda = lambda, stdToo = FALSE)$yt, tfarloc, tfarsca)
        pn <- rep(0, length(qs))
        pn[qsIndnz] <- pnorm(zn)
        return(pn)
      }
    }
    funcsList[[g]] <- tfunc(farloc, farsca, lambda, tfarloc, tfarsca)
    probsMatrix[,g] <- probs
  }
  return(list(lDACs = lDACValues, farness = farnessFinal, gfarness = apply(probsMatrix, 1, min), farscales = farscales, funcsList = funcsList))
}

testlDAC <- function(trainlDACList, distMatrix, trueClass) {
  #Apply functions with parameters computed in trainlDAC
  #Assumes distance matrix has column names
  
  #Returns: 
  #lDAC and farness of test dataset using funcsList and farscales from trainPac
  #gfarness - global farness to find outlier observation to every class
  
  uniqueClass <- sort(unique(trueClass))
  numClass <- length(uniqueClass)
  lDACValues <- c()
  #Compute lDAC values
  for (i in seq(1, nrow(distMatrix))) {
    distTrue <- distMatrix[i, trueClass[i]]
    distPredicted <- min(distMatrix[i, -which(colnames(distMatrix) == trueClass[i])])
    lDACValues <-  c(lDACValues, exp(-distPredicted/2)/(exp(-distPredicted/2) + exp(-distTrue/2)) )
  }
  farness <- rep(0, ncol(distMatrix))
  farscales <- trainlDACList$farscales
  funcsList <- trainlDACList$funcsList
  #Yeo-Johnson transform
  fig <- scale(distMatrix, center = FALSE, scale = farscales)
  for (g in seq(1, numClass)) {
    far2probTemp <- funcsList[[g]]
    fig[, g] <- far2probTemp(fig[, g])
  }
  for (g in seq(1, numClass)) {
    classInd <- which(trueClass == uniqueClass[g])
    farness[classInd] <- fig[classInd, g]
  }
  return(list(lDACs = lDACValues, farness = farness, gfarness = apply(fig, 1, min), fig = fig))
}

computeProbs <- function(trueClass, predictedClass, gfarness = NULL, plotGfarness = FALSE, tau = 0.99) {
  #Compute Probability matrix
  uniqueClass <- sort(unique(trueClass))
  true <- factor(trueClass, levels = uniqueClass)
  predicted <- factor(predictedClass, levels = uniqueClass)
  if (plotGfarness) {
    gfarInd <- which(gfarness > tau)
    if (length(gfarInd) > 0) {
      predictedClass <- as.character(predictedClass)
      predictedClass[gfarInd] <- "outl"
      levName <- c(as.character(uniqueClass), "outl")
      predicted <- factor(predictedClass, levels = levName)
    }
  }
  probMatrix <- table(true, predicted)
  probMatrix <- probMatrix/rowSums(probMatrix)
  return(probMatrix)
}

plotMosaic <- function(trueClass, predictedClass, fillColor, gfarness = NULL, plotGfarness = FALSE, tau = 0.99) {
  #Mosaic plot of classes. fillColor is a vector of colors for the classes.
  #If plotGfarness is TRUE the gfarness observations are shown on top in black "#4F4F4F".
  probMatrix <- computeProbs(trueClass, predictedClass, gfarness, plotGfarness, tau)
  uniqueClass <- sort(unique(trueClass))
  numClass <- length(uniqueClass)
  
  permuteColor <- function(n) {
    #Colors grow from bottom to top
    cmat <- matrix(0, n, n)
    for (i in seq(1, n)) {
      cmat[1, i] <- i
      cmat[-1, i] <- seq(1, n)[-i]
    }
    return(cmat)
  }
  
  pmatColor <- permuteColor(numClass)
  
  if (plotGfarness && !is.null(gfarness)) {
    gfarInd <- which(gfarness > tau)
    if (length(gfarInd) > 0) {
      fillColor <- c(fillColor, "#4F4F4F")
      pmatColor <- rbind(pmatColor, rbind(rep(numClass + 1, numClass)))  
    } else {
      plotGfarness <- FALSE
    }
  }
  
  counts <- colorsList <- list()
  for (i in seq(1, numClass)) {
    colorAux <- countAux <- c()
    for (j in seq(1, numClass + plotGfarness)) {
      count <- round(probMatrix[i,j]*100)
      countAux <- c(countAux, count)
      colorAux <- c(colorAux, fillColor[j])
    }
    countInd <- which(countAux != 0)
    pmatColorAux <- pmatColor[,i][pmatColor[,i] %in% countInd]
    counts[[i]] <- countAux[pmatColorAux]
    colorsList[[i]] <- colorAux[pmatColorAux]
  }
  
  factorTrueClass <- factor(trueClass, levels = uniqueClass)
  widths <- table(factorTrueClass)/length(trueClass)
  pos <- 0.5 * (cumsum(widths) + cumsum(c(0, widths[-length(widths)])))
  
  gplot <- ggplot2::ggplot()
  
  for (i in seq(1, numClass)) {
    posi <- pos[i]
    counti <- counts[[i]]
    filli <- colorsList[[i]]
    gplot <- gplot +
      geom_bar(aes_(x = posi, y = counti), 
               fill = filli, stat = "identity", 
               width = widths[i], color = "white", 
               position = position_fill())
  }
  
  gplot <- gplot + 
    xlab("True Class") + ylab("Predicted Class") +
    scale_x_continuous(breaks = pos,
                       labels = uniqueClass,
                       expand = c(0.02, 0)) +
    scale_y_continuous(breaks = c(0,1),
                       expand = c(0.02, 0)) +
    theme(axis.text.x = element_text(size = 14, face = "italic"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 16, face = "bold"), 
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.y = element_text(size = 14), 
          axis.title.y = element_text(size = 17, face = "bold"))
  
  plot(gplot)
}

plotClassmap <- function(trueClass, predictedClass, lDACs, farness, fillColor, pointSize, gfarness, tau = 0.99) {
  #Plots of lDAC and farness per class
  #Assume the lDACs and farness have been calculated previously
  #The gfarness observations are shown using a triangle.
  uniqueClass <- sort(unique(trueClass))
  nTrue <- length(uniqueClass)
  nlDAC <- length(lDACs)
  groupInfo <- data.frame(matrix(0, ncol = 4, nrow = nlDAC))
  groupInfo[,1] <- lDACs
  groupInfo[,2] <- farness
  groupInfo[,3] <- trueClass
  groupInfo[,4] <- predictedClass  
  
  colornames <- setNames( fillColor, uniqueClass)
  
  scale_new <- c(0, 0.5, 0.7, 0.9, 0.99, 0.999, 1)
  scale_new_length <- c(0.5, 0.7 - 0.5, 0.9 - 0.7, 0.99 - 0.9, 0.999 - 0.99, 1 - 0.999)
  scale_real <- c(0, 2/9, 1/3, 4/9, 2/3, 7/9)
  scale_real_length <- c(2/9, 1/9, 1/9, 2/9, 1/9, 2/9)
  
  #Scale tau
  tau_index <- which((tau <= scale_new) == 1)[1] - 1
  tau_transform <- (tau - scale_new[tau_index])/scale_new_length[tau_index] * scale_real_length[tau_index] + scale_real[tau_index]
  
  for (cl in seq(1, nTrue)) {
    classInd <- which(trueClass == uniqueClass[cl])
    predClass <- predictedClass[classInd]
    uniquePred <- sort(unique(predClass))
    predAux <- c()
    for (pInd in seq(1, length(uniquePred))) {
      predAux <- c(predAux, which(uniqueClass == uniquePred[pInd]))
    }
    uniqueColor <- fillColor[predAux]
    classFarness <- farness[classInd]
    classGfarness <- (gfarness[classInd] > tau)
    classlDACs <- lDACs[classInd]
    pointShape <- rep(16, length(classInd))
    pointShape[classGfarness] <- 17
    colorVec <- c()
    for (i in classInd) {
      colorVec <- c(colorVec, fillColor[which(uniqueClass == predictedClass[i])])
    }
    classFarnessTransform <- c()
    for (j in seq(1, length(classFarness))) {
      classFarnessj <- classFarness[j]
      for (iter in seq(1, 6)) {
        if (classFarnessj <= scale_new[iter + 1]) {
          classFarnessj <- (classFarnessj - scale_new[iter])/scale_new_length[iter] * scale_real_length[iter] + scale_real[iter]
          break
        }
      }
      classFarnessTransform <- c(classFarnessTransform, classFarnessj)
    }
    
    gplot <- ggplot(data = data.frame(classlDACs = classlDACs, classFarnessTransform = classFarnessTransform)) +
      geom_point(aes(x = classFarnessTransform, y = classlDACs, color = colorVec, size = "a"), shape = pointShape) +
      geom_hline(yintercept = 1/2) + 
      geom_vline(aes(xintercept = tau_transform, linetype = "dashed")) + 
      scale_linetype_manual(values = "dashed", guide = "none") +
      scale_size_manual(values = pointSize, guide = "none") +
      scale_color_manual(name = "Class", labels = uniquePred,
                         values = eval(parse(text = paste("c(", toString(paste("\'" , uniqueColor, "\'", "=" , "\'", uniqueColor,"\'", sep = "")), ")")) )
                         ) + 
      xlab(bquote("Farness (" * tau==.(tau) * ")")) +
      ylab(expression("\u2113" * "DAC")) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
      theme(axis.title.y = element_text(size = 15, face = "bold"), 
            axis.title.x = element_text(size = 15, face = "bold")) + 
      scale_x_continuous(breaks = c(0, 2/9, 4/9, 2/3, 7/9, 1), 
                         labels = c("0", "0.5", "0.9", "0.99", "0.999", "1"), 
                         minor_breaks = seq(0, 1, 1/9))
    
    print(gplot)
  }
  colnames(groupInfo) <- c("lDAC", "Farness", "True", "Predicted")
  return(groupInfo)
}

plotSilhouette <- function(trueClass, lDACs, plot = TRUE, fillColor) {
  # Determine assigned class: smallest distance = predicted class
  assignedLabels <- trueClass
  ord <- order(assignedLabels, decreasing = TRUE)
  
  # Build silhouette-like data frame (obs in original 1:n order but clustered by 'ord')
  dataSilh <- data.frame(
    silhWidth = 1 - 2 * lDACs[ord],
    cluster = factor(assignedLabels[ord]),
    obs = factor(1:length(assignedLabels))
  )
  
  dataSilhSorted <- dataSilh
  clusterLevels <- levels(dataSilh$cluster)
  for (cl in clusterLevels) {
    idx <- which(dataSilh$cluster == cl)
    if (length(idx) > 1) {
      # sort ascending and assign to the indices in their original order
      dataSilhSorted$silhWidth[idx] <- sort(dataSilh$silhWidth[idx], decreasing = FALSE)
    } else {
      dataSilhSorted$silhWidth[idx] <- dataSilh$silhWidth[idx]
    }
  }
  
  dataSilh <- dataSilhSorted
  
  # Compute average silhouette per class
  silhSummary <- dataSilh %>%
    group_by(cluster) %>%
    summarise(silhMean = mean(silhWidth))
  
  # Overall average
  overallMean <- mean(dataSilh$silhWidth)
  
  if (plot) {
    # Uniform colors per class
    nClass <- nlevels(dataSilh$cluster)
    
    # Create vertical silhouette plot
    gplot <- ggplot(dataSilh, aes(y = obs, x = silhWidth, fill = cluster)) +
             geom_col(width = 1) +
             geom_vline(xintercept = 0, color = "gray40", linetype = "dashed") +
             scale_fill_manual(name = "Class", values = fillColor, labels = sort(unique(trueClass))) +
             labs(title = "", x = "Silhouette width s(j)", y = NULL) +
             theme_minimal(base_size = 13) +
             theme(
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "plain", size = 13),
               axis.title.x = element_text(size = 13, , face = "plain", margin = margin(t = 15)),
               plot.margin = margin(10, 20, 60, 20)
             )
    
    # Prepare mean silhouette labels
    silhSummary$label <- paste0("s\u0304(", silhSummary$cluster, ") = ",
                                sprintf("%.2f", silhSummary$silhMean))
    
    # Position for each class label (top of block)
    labelPos <- dataSilh %>%
      group_by(cluster) %>%
      summarise(yPos = max(as.numeric(obs)))
    
    silhSummary <- left_join(silhSummary, labelPos, by = "cluster")
    
    # Place mean silhouette labels slightly LEFT of x = 0
    gplot <- gplot + geom_text(
             data = silhSummary,
             aes(x = -0.05, y = yPos, label = label),
             color = "black", hjust = 1,
             size = 4.2, fontface = "plain"
             )
    
    print(gplot)
    grid.text(sprintf("Overall average silhouette width: %.2f", overallMean), x = unit(0.45, "npc"), y = unit(.1, "npc"), just = "centre")
  }
  
  return(list( dataSilh = dataSilh, silhSummary = silhSummary, overallMean = overallMean))
}

plotFarness <- function(trueClass, predictedClass, farness, fillColor, pointSize, gfarness, tau = 0.95) {
  # Plot farness vs observation index
  # Outliers (gfarness > tau) marked as triangles
  # Vertical lines separate true classes (placed between groups)
  
  n <- length(trueClass)
  obsIndex <- 1:n
  
  # Shapes: circle for normal, triangle for outliers
  pointShape <- ifelse(gfarness > tau, 17, 16)
  
  # Colors depend on predicted class
  uniqueClass <- sort(unique(trueClass))
  nClass <- length(uniqueClass)
  predClassFactor <- factor(predictedClass, levels = uniqueClass)
  
  df <- data.frame(
    obsIndex = obsIndex,
    farness = farness,
    predClassFactor = predClassFactor,
    pointShape = factor(pointShape, levels = c(16, 17),
                        labels = c("Regular", "Global Outlier"))
  )
  
  # Compute class boundaries *between* groups
  classBoundaries <- cumsum(table(trueClass))
  vlines <- classBoundaries[-length(classBoundaries)] + 0.5
  
  # Build y breaks
  ybreaks <- sort(unique(c(pretty(range(farness)), tau)))
  
  gplot <- ggplot(df, aes(x = obsIndex, y = farness,
                          color = predClassFactor, shape = pointShape)) +
    geom_point(size = pointSize) +
    # Horizontal tau line
    geom_hline(yintercept = tau, linetype = "solid", color = "black",
               size = 1) +
    # Vertical separators
    geom_vline(xintercept = vlines, linetype = "dotted", color = "black", 
               size = 1) +
    scale_shape_manual(name = "Observation type",
                       values = c("Regular" = 16, "Global Outlier" = 17)) +
    scale_color_manual(name = "Predicted class",
                       values = fillColor) +
    scale_y_continuous(
      name = "Farness",
      breaks = ybreaks,
      labels = function(x) {
        ifelse(abs(x - tau) < .Machine$double.eps^0.5,
               expression(bold(tau)), x)
      }
    ) +
    xlab("Observation index") +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      axis.line = element_blank(),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  print(gplot)
}
