source("IFDA.R")
source("lDAC.R")

#Load Credit Card Dataset
load("credicard_CR.RDATA")

creditDims <- dim(CreditCard_CR[,1:5])
creditTrueClass <- c(rep("User 1", 12), rep("User 2", 12), rep("User 3", 12))
deltaT <- 1/24
credVarib.C <- c(1,4)

creditSymbDataset <- dataSplit(CreditCard_CR[,credVarib.C], CreditCard_CR[,(credVarib.C+5)], creditTrueClass, 1, seed = 3)

#Apply Ifda
creditSymbIfda <- ifda(deltaT, creditSymbDataset$trainC, creditSymbDataset$trainR, 
                       creditSymbDataset$trainC, creditSymbDataset$trainR, 
                       creditSymbDataset$trainClass, creditSymbDataset$trainClass, 
                       restriction = "uncorrelated", s=3)

##Confusion matrix
table(creditSymbDataset$trainClass, creditSymbIfda$predictedClass)
sum(diag(table(creditSymbDataset$trainClass, creditSymbIfda$predictedClass)))/length(creditSymbDataset$trainClass)

#Compute distance matrix of train dataset
distMatrixTrainSymb <- mallowsPredict(creditSymbIfda$alphaMatrix, deltaT, 
                                      creditSymbDataset$trainC, creditSymbDataset$trainR, 
                                      creditSymbIfda$muClassC,  creditSymbIfda$muClassR, 
                                      unique(creditSymbDataset$trainClass))$distMatrix

#Train dataset lDAC
creditlDACTrainSymb <- trainlDAC(distMatrixTrainSymb, creditSymbDataset$trainClass)

#Test dataset lDAC
creditlDACTestSymb <- testlDAC(creditlDACTrainSymb, creditSymbIfda$distMatrix, creditSymbDataset$trainClass)

###Plots

#Mosaic Plot (no outliers)
creditSymbPlotMosaic <- plotMosaic(creditSymbDataset$trainClass, creditSymbIfda$predictedClass, 
                                   c("#ff8c1a","#558C58", "#4A90D5"), 
                                   gfarness = creditlDACTestSymb$gfarness,
                                   plotGfarness = FALSE)

#Mosaic Plot (global outliers tau = 0.95)
creditSymbPlotMosaicOut <- plotMosaic(creditSymbDataset$trainClass, creditSymbIfda$predictedClass, 
                                   c("#ff8c1a", "#558C58", "#4A90D5"), 
                                   gfarness = creditlDACTestSymb$gfarness,
                                   plotGfarness = TRUE, tau = 0.95)

#Classmap
creditSymbPlotClassmap <- plotClassmap(creditSymbDataset$trainClass, creditSymbIfda$predictedClass, 
                                       creditlDACTestSymb$lDACs, creditlDACTestSymb$farness, 
                                       c("#ff8c1a", "#558C58", "#4A90D5"), 
                                       3, creditlDACTestSymb$gfarness, tau = 0.95)

#Farness plot
creditSymbPlotFarness <- plotFarness(creditSymbDataset$trainClass,
                                     creditSymbIfda$predictedClass,
                                     creditlDACTestSymb$farness,
                                     c("#ff8c1a","#558C58", "#4A90D5"),
                                     3,
                                     creditlDACTestSymb$gfarness, tau=0.95)

#Silhouette plot
creditSymbPlotSilhouette <- plotSilhouette(trueClass = creditTrueClass,
                                           lDACs = creditlDACTestSymb$lDACs,
                                           plot = TRUE, c("#ff8c1a","#558C58", "#4A90D5"))
