# Learn the model parameter from the synthetic data 
#
# Author: Jun Yu
# Version: Jan 18, 2012
##################################################################

rm(list=ls())

#setwd("C:/Jun_home/workspace/eBird.ODFP")
#setwd("/Users/yujunnokia/Documents/workspace/eBird.ODFP")
setwd("/Users/yujunnokia/workspace/eBird.ODFP")
source("ODFP.R")
source("ODFP.synthData.R")

library("lattice")
library("Matrix")
library("glmnet")

######################
# experiment settings
######################
nTrSites <- 500  # number of training sites
nTeSites <- 500  # number of testing sites
nVisits <- 3  # number of visits to each site
nOccCovs <- 5  # number of occupancy covariates
nDetCovs <- 5  # number of detection covariates
nParams <- nOccCovs + nDetCovs * 2 + 3 # total number of paramters. 2 is the two intercept terms from occupancy and detection parts
nRandomRestarts <- 1  # number of random restarts  

#################
# regularization
#################
regType <- 2 # regularization types: 0 for none, 1 for L1, 2 for L2
lambda <- lambdaO <- lambdaD <- 0.01  # regularization paramters

#######################
# Set model parameters
#######################
alpha <- c(1,rnorm(nOccCovs)*3)
gamma <- c(1,rnorm(nDetCovs, mean = 1, sd = 1)*3)
eta <- c(1,rnorm(nDetCovs, mean = -1, sd = 1)*3)
covs <- rnorm(1000*nDetCovs, mean = 0.5, sd = 1)
dim(covs) <- c(1000,nDetCovs)
covs <- cbind(1,covs)
cat("mean true detection prob is",mean(Logistic(covs %*% gamma)),"\n",sep=" ")
cat("mean false detection prob is",mean(Logistic(covs %*% eta)),"\n",sep=" ")


########################
# generate testing data
########################
teVisits <- array(nVisits, nTeSites)
teData <- GenerateData(nTeSites,teVisits,alpha,gamma,eta)
teDetHists <- teData$detHists
teOccCovs <- teData$occCovs
teDetCovs <- teData$detCovs
teTrueOccs <- teData$trueOccs

#########################
# generate training data
#########################
trVisits <- array(0, c(nTrSites,1))
for (i in 1:nTrSites) {
    isMultipleVisits <- runif(1) < 0.5
    if (isMultipleVisits == TRUE) {
        trVisits[i] <- round(runif(1, min=2, max=nVisits))
    } else {
        trVisits[i] <- 1
    }
} # i
trData <- GenerateData(nTrSites,trVisits,alpha,gamma,eta)
trDetHists <- trData$detHists
trOccCovs <- trData$occCovs
trDetCovs <- trData$detCovs
trTrueOccs <- trData$trueOccs

#####################
# get Bayes rates
#####################
{
    # get occupancy rate and detection rate
    teOccProb <- array(0, c(nTeSites, 1))
    teTrueDetProb <- array(0, c(nTeSites, nVisits))
    teFalseDetProb <- array(0, c(nTeSites, nVisits))
    predDetHists <- array(0, c(nTeSites, nVisits))
    for (i in 1:nTeSites) {
        teOccProb[i] <- Logistic(teOccCovs[i,] %*% alpha)
        for (t in 1:teVisits[i]) {
            teTrueDetProb[i,t]  <- Logistic(teDetCovs[i,t,] %*% gamma)
            teFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% eta)
             
            if (round(teOccProb[i]) == 1) {
                predDetHists[i,t] <- round(teTrueDetProb[i,t])
            } else {
                predDetHists[i,t] <- round(teFalseDetProb[i,t])				
            }
        } # t
    } # i
    bayesOcc <- sum(round(teOccProb) == teTrueOccs) / nTeSites
    bayesDet <- sum(sum(predDetHists == teDetHists)) / (sum(teVisits))
    cat("bayes occupancy rate is ", bayesOcc, "\n")
    cat("bayes detection rate is ", bayesDet, "\n")
}

#############
# ODFP model
#############
{
    # run ODFP
    params <- RandomRestartEM(trDetHists,trOccCovs,trDetCovs,trVisits,regType,lambdaO,lambdaD,nRandomRestarts)
    alphaODFP <- params$alpha
    gammaODFP <- params$gamma
    etaODFP <- params$eta
    
    # compute occupancy rate and detection rate
    teOccProb <- array(0,c(nTeSites,1))
    teTrueDetProb <- array(0,c(nTeSites,nVisits))
    teFalseDetProb <- array(0,c(nTeSites,nVisits))
    predDetHists <- array(0,c(nTeSites,nVisits))
    for (i in 1:nTeSites) {
        teOccProb[i] <- Logistic(teOccCovs[i,] %*% alphaODFP)
        for (t in 1:teVisits[i]) {
            teTrueDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% gammaODFP)
            teFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% etaODFP)
            
            if (round(teOccProb[i]) == 1) {
                predDetHists[i,t] <- round(teTrueDetProb[i,t])
            } else {
                predDetHists[i,t] <- round(teFalseDetProb[i,t])				
            }
        } # t
    } # i
    modelOcc <- sum(round(teOccProb) == teTrueOccs) / nTeSites
    modelDet <- sum(sum(predDetHists == teDetHists)) / (sum(teVisits))
    cat("------------------------------\n")
    cat("bayes occupancy rate is ",bayesOcc,"\n")
    cat("model occupancy rate is ",modelOcc,"\n")
    cat("bayes detection rate is ",bayesDet,"\n")
    cat("model detection rate is ",modelDet,"\n")
    
    
    # predict Z on test data
    teBayesOccProb <- array(0,c(nTeSites,1))
    teModelOccProb <- array(0,c(nTeSites,1))
    for (i in 1:nTeSites) {
        teBayesOccProb[i] <- PredictOcc(c(alpha,gamma,eta),teOccCovs[i,],teDetCovs[i,,],teDetHists[i,],teVisits[i]) 
        teModelOccProb[i] <- PredictOcc(c(alphaODFP,gammaODFP,etaODFP),teOccCovs[i,],teDetCovs[i,,],teDetHists[i,],teVisits[i])    
    } # i
    trueOcc = sum(round(teBayesOccProb) == teTrueOccs) / nTeSites
    predOcc = sum(round(teModelOccProb) == teTrueOccs) / nTeSites
    cat("------------------------------\n")
    cat("True occupancy prediction is ",trueOcc,"\n")
    cat("Model occupancy prediction is ",predOcc,"\n")
    
    # predict Y on test data
    trueDetHists <- array(0,c(nTeSites,nVisits))
    modelDetHists <- array(0,c(nTeSites,nVisits))
    for (i in 1:nTeSites) {
        for (t in 1:teVisits[i]) {
            trueDetHists[i,t] <- PredictDet(c(alpha,gamma,eta),teOccCovs[i,],teDetCovs[i,t,]) 
            modelDetHists[i,t] <- PredictDet(c(alphaODFP,gammaODFP,etaODFP),teOccCovs[i,],teDetCovs[i,t,]) 
        }
    }
    trueDet <- sum(sum(round(trueDetHists) == teDetHists)) / (sum(teVisits))
    predDet <- sum(sum(round(modelDetHists) == teDetHists)) / (sum(teVisits))
    cat("True detection prediction is ",trueDet,"\n")
    cat("Model detection prediction is ",predDet,"\n")
    
    # compute MSE
    MSE <- sum(sum((alpha-alphaODFP)^2) + sum((gamma-gammaODFP)^2) + sum((eta-etaODFP)^2)) / nParams
    cat("------------------------------\n")
    cat("MSE is ",MSE,"\n")
}

############
## OD model
############
#{
#    source("../eBird.OD/OD.R")
#    params <- RandomRestartEM(trDetHists,trOccCovs,trDetCovs,trVisits,regType,lambdaO,lambdaD,nRandomRestarts)
#    alphaOD <- params$alpha
#    betaOD <- params$beta
#    
#    cat("------------------------------\n")
#    
#    # predict Z on test data
#    teModelOccProb <- array(0,c(nTeSites,1))
#    for (i in 1:nTeSites) {
#        teModelOccProb[i] <- PredictOcc(c(alphaOD,betaOD),teOccCovs[i,],teDetCovs[i,,],teDetHists[i,],teVisits[i]) 
#    } # i
#     predOcc <- sum(round(teModelOccProb) == teTrueOccs) / nTeSites
#    cat("Model occupancy prediction is ",predOcc,"\n")
#    
#    # predict Y on test data
#    teModelDetHists <- array(0,c(nTeSites,nVisits))
#    for (i in 1:nTeSites) {
#        for (t in 1:teVisits[i]) {
#            teModelDetHists[i,t] <- PredictDet(c(alphaOD,betaOD),teOccCovs[i,],teDetCovs[i,t,]) 
#        }
#    }
#    predDet <- sum(sum(round(teModelDetHists) == teDetHists)) / (sum(teVisits))
#    cat("Model detection prediction is ",predDet,"\n")
#}
#
#
##################################################
## run logistic regression to predict detection Y
##################################################
#{
#    trOccDetCovs <- NULL 
#    trDetections <- NULL 
#    for (i in 1:nTrSites) {
#        for (t in 1:trVisits[i]) {
#            trOccDetCovs <- rbind(trOccDetCovs,c(trOccCovs[i,2:(nOccCovs+1)],trDetCovs[i,t,2:(nDetCovs+1)]))
#            trDetections <- rbind(trDetections,trDetHists[i,t])
#        }
#    }
#    
#    # train GLMNET with L2 norm
#    LRM <- glmnet(trOccDetCovs,trDetections,family="binomial",alpha=0)
#    
#    # get occupancy rate
#    teOccDetCovs <- NULL 
#    teDetections <- NULL 
#    for (i in 1:nTeSites) {
#        for (t in 1:teVisits[i]) {
#            teOccDetCovs <- rbind(teOccDetCovs,c(teOccCovs[i,2:(nOccCovs+1)],teDetCovs[i,t,2:(nDetCovs+1)]))
#            teDetections <- rbind(teDetections,teDetHists[i,t])
#        }
#    }
#    
#    # predict
#    predDetHists <- predict(LRM,type="response",newx=teOccDetCovs,s=lambda)
#    predDetHists[predDetHists >= 0.5] <- 1
#    predDetHists[predDetHists < 0.5] <- 0
#    
#    LRMDet <- sum(predDetHists == teDetections) / sum(teVisits)
#    cat("LRM detection rate is ",LRMDet,"\n")
#}
