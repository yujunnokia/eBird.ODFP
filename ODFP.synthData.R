# TODO: generate synthetic data from the ODFP model
#
# Author: Jun Yu
# Version: Jan 18, 2012
###############################################################################

source("ODFP.R")

#
# generate synthetic data from ODFP model
#
# Args:
#   nSites: number of sites
#	visits: a vector recording number of visits per site
#	alpha: occupancy parameter
#	gamma: detection parameter for true detection
#	eta: detection parameter for false detection
#
# Returns:
#   synthetic data generated from the model including detection history, occupancy 
#	covaraites, detection covariates and true occupancy status.
#
GenerateData <- function(nSites,visits,alpha,gamma,eta) 
{
    nOccCovs <- length(alpha)
    nDetCovs <- length(gamma)
    nVisits  <- max(visits)
    
    # generate occupancy and detection covariates
    occCovs <- rnorm(nSites * nOccCovs, mean = 0, sd = 1)
    dim(occCovs) <- c(nSites, nOccCovs)
    occCovs[,1] <- array(1, c(nSites, 1))
    
    detCovs <- rnorm(nSites * nVisits * nDetCovs, mean = 0.5, sd = 1)
    dim(detCovs) <- c(nSites, nVisits, nDetCovs)
    detCovs[,,1] <- array(1, c(nSites, nVisits, 1))
    
    trueOccs <- array(0, c(nSites, 1))
    trueOccs <- runif(nSites) < Logistic(occCovs %*% alpha)
    
    detHists <- array(0, c(nSites, nVisits))
    for (i in 1:nSites) {
        for (t in 1:visits[i]) {
            # if the bird is truely present and trueObs is 1, then detHists is 1. O.w. detHists is 0.
            # if the bird is truely absent and falseObs is 1, then detHists is 1. O.w. detHists is 0.
            if (trueOccs[i] == 1) {
                isTrueDetected <- runif(1) < Logistic(detCovs[i,t,] %*% gamma)
                if (isTrueDetected == 1) {
                    detHists[i,t] <- 1
                } 
            } else {
                isFalseDetected <- runif(1) < Logistic(detCovs[i,t,] %*% eta)
                if (isFalseDetected == 1) {
                    detHists[i,t] <- 1
                } 
            } 
        } # t
    } # i
    
    retData <- list(detHists=detHists, occCovs=occCovs, detCovs=detCovs, trueOccs=trueOccs)
    return(retData)
}

