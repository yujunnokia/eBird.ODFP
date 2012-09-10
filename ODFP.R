# TODO: implementation of ODFP model
#
# Author: Jun Yu
# Version: Jan 18, 2012
###############################################################################


#####################
# Utility functions #
#####################

#
# Compute logistic function
#
# Args:
#   x: input to logistic function
#
# Returns:
#   Output of logistic function
#
Logistic <- function(x) 
{
    y <- 1 / ( 1 + exp(-x))    
    return(y)
}


#
# Implement log(exp(a) + exp(b))
#
# Args:
#   a: first input argument 
#	b: second input argument
#
# Returns:
#   output log(exp(a) + exp(b))
#
LogSumExp <- function(a,b) 
{
    c <- -Inf
    if (b < a) {
        c <- a + log(1 + exp(b - a))
    } else {
        if (a == -Inf && b == -Inf) {
            c <- -Inf
        } else {
            c <- b + log(1 + exp(a - b))
        }
    }    
    return(c)
}


#
# Implement log(exp(a) - exp(b))
#
# Args:
#   a: first input argument 
#	b: second input argument
#
# Returns:
#   output log(exp(a) - exp(b))
#
LogDiffExp <- function(a,b) 
{
    c <- -Inf    
    if (b < a) {
        c <- a + log(1 - exp(b - a))
    } else if (a == b) {
        c <- -Inf
    } else {
        warning("LogDiffExp output -inf.\n")
    }  
    return(c)
}


#########################
# OD.FP model functions #
#########################

#
# E-step: estimate the probability of Z given the model parameters
#
# Args:
#   params: a vector of occupancy and detection parameters
#	Y: observation matrix of size nSites * nVisits 
#	Xo: occupancy covariate matrix of size nSite * nOccCovs
#	Xd: detection covariate matrix of size nSite * nVisits * nDetCovs
#	visits: visit vecter recording number of visits to each site
#
# Returns:
#   probability of expected occupancy at each site
#
ExpectationStep <- function(params,Y,Xo,Xd,visits) 
{
    nSites   <- dim(Y)[1]  # number of sites
    nVisits  <- dim(Y)[2]  # number of visits
    nOccCovs <- dim(Xo)[2]  # number of occupancy covs
    nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    alpha <- params[1:nOccCovs]  # get alpha
    gamma <- params[(nOccCovs + 1):(nOccCovs + nDetCovs)]  # get gamma
    eta   <- params[(nOccCovs + nDetCovs + 1):(nOccCovs + 2 * nDetCovs)]  # get eta
    
    probOcc <- Logistic(Xo %*% alpha)
    probTrueDet <- array(0, c(nSites, nVisits))
    probFalseDet <- array(0, c(nSites, nVisits))
    for (t in 1:nVisits) {
        probTrueDet[,t] <- Logistic(Xd[,t,] %*% gamma)
        probFalseDet[,t] <- Logistic(Xd[,t,] %*% eta)
    } # t
    logProbOcc <- log(probOcc)
    logProbTrueDet <- log(probTrueDet)
    logProbFalseDet <- log(probFalseDet)
    
    # compute probability of expected occupancy
    probExpectedOcc <- rep(0,nSites)
    for (i in 1:nSites) {
        if (logProbOcc[i] == -Inf) {
            probExpectedOcc[i] <- 0
            next
        }
        if (LogDiffExp(0,logProbOcc[i]) == -Inf) {
            probExpectedOcc[i] <- 1
            next
        }
        
        temp1 <- 0
        temp2 <- 0
        for (t in 1:visits[i]) {
            if (Y[i,t] != 0 || logProbTrueDet[i,t] != -Inf) {
                temp1 <- temp1 + Y[i,t] * logProbTrueDet[i,t] 
            } 
            if ((1 - Y[i,t]) != 0 || LogDiffExp(0, logProbTrueDet[i,t]) != -Inf) {
                temp1 <- temp1 + (1 - Y[i,t]) * LogDiffExp(0, logProbTrueDet[i,t])
            }
            
            if (Y[i,t] != 0 || logProbFalseDet[i,t] != -Inf) {
                temp2 <- temp2 + Y[i,t] * logProbFalseDet[i,t] 
            } 
            if ((1 - Y[i,t]) != 0 || LogDiffExp(0, logProbFalseDet[i,t]) != -Inf) {
                temp2 <- temp2 + (1 - Y[i,t]) * LogDiffExp(0, logProbFalseDet[i,t])
            }
            
#			temp1 <- temp1 + Y[i,t]*logProbTrueDet[i,t] + (1-Y[i,t])*LogDiffExp(0,logProbTrueDet[i,t])
#			temp2 <- temp2 + Y[i,t]*logProbFalseDet[i,t] + (1-Y[i,t])*LogDiffExp(0,logProbFalseDet[i,t])
        } # t
        probExpectedOcc[i] <- exp(logProbOcc[i] + temp1) / (exp(logProbOcc[i] + temp1) + exp(LogDiffExp(0,logProbOcc[i]) + temp2))
        
        if (is.na(probExpectedOcc[i])) {
            stop(paste("probExpectedOcc of site", i, "is NaN", sep=" "))
        }
    } # i
    
    return(probExpectedOcc)
}

#
# compute expected joint log-likelihood (EJLL)
#
# Args:
#   params: a vector of occupancy and detection parameters
#	Y: observation matrix of size nSites * nVisits 
#	Xo: occupancy covariate matrix of size nSite * nOccCovs
#	Xd: detection covariate matrix of size nSite * nVisits * nDetCovs
#	visits: visit vecter recording number of visits to each site
#	probExpectedOcc: probablity of expected occupancy Z_i
#	regType: regularization type
#	lambdaO: regularization parameter for occupancy component 
#	lambdaD: regularization parameter for detection component
#
# Returns:
#   expected joint log likelihood
#
ComputeEJLL <- function(params,Y,Xo,Xd,visits,probExpectedOcc,regType,lambdaO,lambdaD) 
{
    nSites   <- dim(Y)[1]  # number of sites
    nVisits  <- dim(Y)[2]  # number of visits
    nOccCovs <- dim(Xo)[2]  # number of occupancy covs
    nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    alpha <- params[1:nOccCovs]  # get alpha
    gamma <- params[(nOccCovs + 1):(nOccCovs + nDetCovs)]  # get gamma
    eta   <- params[(nOccCovs + nDetCovs + 1):(nOccCovs + 2 * nDetCovs)]  # get eta
    
    probOcc <- Logistic(Xo %*% alpha)
    probTrueDet <- array(0, c(nSites, nVisits))
    probFalseDet <- array(0, c(nSites, nVisits))
    for (t in 1:nVisits) {
        probTrueDet[,t] <- Logistic(Xd[,t,] %*% gamma)
        probFalseDet[,t] <- Logistic(Xd[,t,] %*% eta)
    } # t
    logProbOcc <- log(probOcc)
    logProbTrueDet <- log(probTrueDet)
    logProbFalseDet <- log(probFalseDet)
    
    eJLL <- 0
    for (i in 1:nSites) {
        temp1 <- 0
        temp2 <- 0
        for (t in 1:visits[i]) {
            if (Y[i,t] != 0 || logProbTrueDet[i,t] != -Inf) {
                temp1 <- temp1 + Y[i,t] * logProbTrueDet[i,t] 
            } 
            if ((1 - Y[i,t]) != 0 || LogDiffExp(0, logProbTrueDet[i,t]) != -Inf) {
                temp1 <- temp1 + (1 - Y[i,t]) * LogDiffExp(0, logProbTrueDet[i,t])
            }
            
            if (Y[i,t] != 0 || logProbFalseDet[i,t] != -Inf) {
                temp2 <- temp2 + Y[i,t] * logProbFalseDet[i,t] 
            } 
            if ((1 - Y[i,t]) != 0 || LogDiffExp(0, logProbFalseDet[i,t]) != -Inf) {
                temp2 <- temp2 + (1 - Y[i,t])*LogDiffExp(0, logProbFalseDet[i,t])
            }            
        } # t
        
        if (probExpectedOcc[i] != 0 || (logProbOcc[i] + temp1) != -Inf) {
            eJLL <- eJLL + probExpectedOcc[i] * (logProbOcc[i] + temp1)
        }
        
        # if 1-probExpectedOcc[i] is zero, then it won't contribute to the eJLL even if temp2 is -inf
        if ((1 - probExpectedOcc[i]) != 0 || (LogDiffExp(0, logProbOcc[i]) + temp2) != -Inf) {
            eJLL <- eJLL + (1 - probExpectedOcc[i]) * (LogDiffExp(0, logProbOcc[i]) + temp2)
        }
    } # i
    
    if (regType == 1) {
        eJLL <- eJLL - lambdaO*sum(abs(alpha[2:length(alpha)])) - lambdaD*sum(abs(gamma[2:(length(gamma))])) - lambdaD*sum(abs(eta[2:(length(eta))]))
    }
    if (regType == 2) {
        eJLL <- eJLL - 0.5*lambdaO*sum(alpha[2:length(alpha)]^2) - 0.5*lambdaD*sum(gamma[2:(length(gamma))]^2) - 0.5*lambdaD*sum(eta[2:(length(eta))]^2)
    }
    
    negEJLL <- -eJLL
    if (is.na(eJLL)) {
        stop("eJLL is na...\n")
    }
    return(negEJLL)
}

#
# compute derivatives of parameter w.r.t. EJLL
#
# Args:
#   params: a vector of occupancy and detection parameters
#	Y: observation matrix of size nSites * nVisits 
#	Xo: occupancy covariate matrix of size nSite * nOccCovs
#	Xd: detection covariate matrix of size nSite * nVisits * nDetCovs
#	visits: visit vecter recording number of visits to each site
#	probExpectedOcc: probablity of expected occupancy Z_i
#	regType: regularization type
#	lambdaO: regularization parameter for occupancy component 
#	lambdaD: regularization parameter for detection component
#
# Returns:
#   Derivatives of parameters
#
ComputeDerivsOfEJLL <- function(params,Y,Xo,Xd,visits,probExpectedOcc,regType,lambdaO,lambdaD) 
{
    nSites   <- dim(Y)[1]  # number of sites
    nVisits  <- dim(Y)[2]  # number of visits
    nOccCovs <- dim(Xo)[2]  # number of occupancy covs
    nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    alpha <- params[1:nOccCovs]  # get alpha
    gamma <- params[(nOccCovs + 1):(nOccCovs + nDetCovs)]  # get gamma
    eta   <- params[(nOccCovs + nDetCovs + 1):(nOccCovs + 2 * nDetCovs)]  # get eta
    
    probOcc <- Logistic(Xo %*% alpha)
    probTrueDet <- array(0, c(nSites, nVisits))
    probFalseDet <- array(0, c(nSites, nVisits))
    for (t in 1:nVisits) {
        probTrueDet[,t] <- Logistic(Xd[,t,] %*% gamma)
        probFalseDet[,t] <- Logistic(Xd[,t,] %*% eta)
    } # t
    
    dQda <- array(0,c(nOccCovs,1))
    dQdg <- array(0,c(nDetCovs,1))
    dQde <- array(0,c(nDetCovs,1))
    for (i in 1:nSites) {
        dQda <- dQda + (probExpectedOcc[i] - probOcc[i]) * Xo[i,]
        
        for (t in 1:visits[i]) {
            dQdg <- dQdg + probExpectedOcc[i] * Xd[i,t,] * (Y[i,t] - probTrueDet[i,t]) 
            dQde <- dQde + (1 - probExpectedOcc[i]) * Xd[i,t,] * (Y[i,t] - probFalseDet[i,t]) 
        } # t
    } # i
    
    if (regType == 1) {
        dQda <- dQda - lambdaO * c(0, sign(alpha[2:length(alpha)]))
        dQdg <- dQdg - lambdaD * c(0, sign(gamma[2:length(gamma)]))
        dQde <- dQde - lambdaD * c(0, sign(eta[2:length(eta)]))
    }
    if (regType == 2) {
        dQda <- dQda - lambdaO * c(0, alpha[2:length(alpha)])
        dQdg <- dQdg - lambdaD * c(0, gamma[2:length(gamma)])
        dQde <- dQde - lambdaD * c(0, eta[2:length(eta)])
    }
    
    derivs <- c(-dQda,-dQdg,-dQde)
    return(derivs)
}

#
# EM
#
# Args:
#	Y: observation matrix of size nSites * nVisits 
#	Xo: occupancy covariate matrix of size nSite * nOccCovs
#	Xd: detection covariate matrix of size nSite * nVisits * nDetCovs
#	visits: visit vecter recording number of visits to each site
#	regType: regularization type
#	lambdaO: regularization parameter for occupancy component 
#	lambdaD: regularization parameter for detection component
#	initialParams: initial parameters
#
# Returns:
#   reestiamted parameters
#
EM <- function(Y,Xo,Xd,visits,regType,lambdaO,lambdaD,initialParams) 
{
    nSites   <- dim(Y)[1]  # number of sites
    nVisits  <- dim(Y)[2]  # number of visits
    nOccCovs <- dim(Xo)[2]  # number of occupancy covs
    nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    params <- initialParams # assign the initial parameters 
    
    # set initial probExpectedOcc
    probExpectedOcc <- array(0.5,nSites)
    for (i in 1:nSites) {
        if (sum(Y[i,]) > 0) {
            probExpectedOcc[i] <- 1
        }
    }
    
    initialEJLL <- -ComputeEJLL(params,Y,Xo,Xd,visits,probExpectedOcc,regType,lambdaO,lambdaD)
    newEJLL <- initialEJLL
    newParams <- params
    diffParams <- 1.0e10
    maxIterations <- 50
    iteration <- 1
    tolerance <- 1.0e-10 #0.01
    
    while (diffParams > tolerance && iteration <= maxIterations) {
        probExpectedOcc <- ExpectationStep(params,Y,Xo,Xd,visits)
        outputs <- optim(params,ComputeEJLL,ComputeDerivsOfEJLL,Y,Xo,Xd,visits,probExpectedOcc,regType,lambdaO,lambdaD,method="BFGS")
        params <- outputs$par
        oldEJLL <- newEJLL
        newEJLL <- -ComputeEJLL(params,Y,Xo,Xd,visits,probExpectedOcc,regType,lambdaO,lambdaD)		

        oldParams <- newParams
        newParams <- params
        diffParams <- sum((newParams-oldParams)^2) / length(newParams)
        
        cat("EM iteration: ", iteration, " EJLL improvement is ", (newEJLL - oldEJLL), 
                "params change is ", diffParams, "\n")
        iteration <- iteration + 1
    }
    
    # check if true detection is higher than false detection
    # if not, reverse the sign of the params
    alpha <- params[1:nOccCovs]  # get alpha
    gamma <- params[(nOccCovs + 1):(nOccCovs + nDetCovs)]  # get gamma
    eta   <- params[(nOccCovs + nDetCovs + 1):(nOccCovs + 2 * nDetCovs)]  # get eta
    
    probTrueDet <- array(0, c(nSites, nVisits))
    probFalseDet <- array(0, c(nSites, nVisits))
    for (t in 1:nVisits) {
        probTrueDet[,t] <- Logistic(Xd[,t,] %*% gamma)
        probFalseDet[,t] <- Logistic(Xd[,t,] %*% eta)
    } # t
    
    if (sum(sum(probTrueDet)) < sum(sum(probFalseDet))) {
        params <- c(-alpha,eta,gamma)
    }
    
    return(params)
}

#
# random restart EM 
#
# Args:
#	Y: observation matrix of size nSites * nVisits 
#	Xo: occupancy covariate matrix of size nSite * nOccCovs
#	Xd: detection covariate matrix of size nSite * nVisits * nDetCovs
#	visits: visit vecter recording number of visits to each site
#	regType: regularization type
#	lambdaO: regularization parameter for occupancy component 
#	lambdaD: regularization parameter for detection component
#	nRandomRestart: number of random restarts
#
# Returns:
#   reestiamted parameters
#
RandomRestartEM <- function(Y,Xo,Xd,visits,regType=2,lambdaO=0.01,lambdaD=0.01,nRandomRestart=2) 
{
	# check the input arguments
	if (dim(Y)[1] != dim(Xo)[1] || dim(Y)[1] != dim(Xd)[1] || dim(Y)[1] != length(visits)) {
		stop("Number of sites does not match.")
	}
	if (!is.element(regType,c(0,1,2))) {
		stop("regularization type is invalid. Need to be 0, 1 or 2")
	}
	
    nSites   <- dim(Y)[1]  # number of sites
    nVisits  <- dim(Y)[2]  # number of visits
    nOccCovs <- dim(Xo)[2]  # number of occupancy covs
    nDetCovs <- dim(Xd)[3]  # number of detection covs
    
    nParams <- nOccCovs + nDetCovs * 2
    
    # in iteration 0, we initialize parameters to be all zeros
    initialParams <- rnorm(nParams)
    params <- EM(Y,Xo,Xd,visits,regType,lambdaO,lambdaD,initialParams)
    probExpectedOccs <- ExpectationStep(params,Y,Xo,Xd,visits)
    alpha <- params[1:nOccCovs]
    gamma <- params[(nOccCovs+1):(nOccCovs+nDetCovs)]
    eta   <- params[(nOccCovs+nDetCovs+1):(nOccCovs+nDetCovs*2)]
    bestEJLL <- -ComputeEJLL(params,Y,Xo,Xd,visits,probExpectedOccs,regType,lambdaO,lambdaD)
    cat("Best EJLL is",bestEJLL,"\n")
    
    # random restart
    for (i in 1:nRandomRestart) {
        cat("============================\n")
        cat("randomRestartEM Iteration ",i,"\n")
        
        done <- FALSE
        while (!done) 
        {
            initialParams <- rnorm(nParams)
            dim(initialParams) <- c(nParams,1)
            
            # run EM
            r <- try(params <- EM(Y,Xo,Xd,visits,regType,lambdaO,lambdaD,initialParams))
            done <- !inherits(r, "try-error")
            cat("done is",done,"\n")
        }
        
        # compute probExpectedOcc
        probExpectedOccs <- ExpectationStep(params,Y,Xo,Xd,visits)
        
        # compute expected joint log likelihood
        newEJLL <- -ComputeEJLL(params,Y,Xo,Xd,visits,probExpectedOccs,regType,lambdaO,lambdaD)
        cat("New ELL is",newEJLL,"\n")
        
        if (newEJLL > bestEJLL) {
            alpha <- params[1:nOccCovs]
            gamma <- params[(nOccCovs+1):(nOccCovs+nDetCovs)]	
            eta   <- params[(nOccCovs+nDetCovs+1):(nOccCovs+nDetCovs*2)]	
            
            bestEJLL <- newEJLL
        }
    }

    print(alpha) 
    print(gamma) 
    print(eta)

    retData <- list(alpha = alpha, gamma = gamma, eta = eta, bestEJLL = bestEJLL)
    return(retData)
}

##############################################

#
# predict the occupancy of a new site given occupancy & detection covariates and detection history
#
# Args:
#	params: parameters of OD model
#	Xo: occupancy covariate vector of size nOccCovs
#	Xd: detection covariate matrix of size nVisits * nDetCovs
#	Y: observation vector of size nVisits 
#	nVisits: number of visits to this site. It is a scalar
#
# Returns:
#   prediction of site occupancy
#
PredictOcc <- function(params,Xo,Xd,Y,nVisits) 
{
    nOccCovs <- length(Xo) # number of occupancy covs
    nDetCovs <- dim(Xd)[2] # number of detection covs
    
	# check input arguments
	if (length(params) != (nOccCovs + nDetCovs*2)) {
		stop("The number of parameters is not consistent with number of occupancy and detetction covariates")
	}
	
    alpha <- params[1:nOccCovs]
    gamma <- params[(nOccCovs+1):(nOccCovs+nDetCovs)]	
    eta   <- params[(nOccCovs+nDetCovs+1):(nOccCovs+nDetCovs*2)]	
    
    if (!is.null(Y)) {
        logProbOcc1 <- log(Logistic(Xo %*% alpha))
        logProbOcc0 <- LogDiffExp(0,logProbOcc1)
        
        for (t in 1:nVisits) {
            if (Y[t] != 0 || log(Logistic(Xd[t,] %*% gamma)) != -Inf) {
                logProbOcc1 <- logProbOcc1 + Y[t] * log(Logistic(Xd[t,] %*% gamma))    
            }
            if ((1 - Y[t]) != 0 || log(1-Logistic(Xd[t,] %*% gamma)) != -Inf) {
                logProbOcc1 <- logProbOcc1 + (1-Y[t]) * log(1-Logistic(Xd[t,] %*% gamma))
            }
             
            if (Y[t] != 0 || log(Logistic(Xd[t,] %*% eta)) != -Inf) {
                logProbOcc0 <- logProbOcc0 + Y[t] * log(Logistic(Xd[t,] %*% eta))    
            }
            if ((1 - Y[t]) != 0 || log(1-Logistic(Xd[t,] %*% eta)) != -Inf) {
                logProbOcc0 <- logProbOcc0 + (1-Y[t]) * log(1-Logistic(Xd[t,] %*% eta))
            }
        }
        
        return(exp(logProbOcc1 - LogSumExp(logProbOcc1,logProbOcc0)))
    } else {
        return(Logistic(Xo %*% alpha))
    }
}


#
# predict the detection of a new visit
#
# Args:
#	params: parameters of OD model
#	Xo: occupancy covariate vector of size nOccCovs
#	Xd: detection covariate vector of size nDetCovs
#
# Returns:
#   prediction of detection
#
PredictDet <- function(params,Xo,Xd) 
{
    nOccCovs <- length(Xo) # number of occupancy covs
    nDetCovs <- length(Xd) # number of detection covs
	
	# check input arguments
	if (length(params) != (nOccCovs + nDetCovs*2)) {
		stop("The number of parameters is not consistent with number of occupancy and detetction covariates")
	}
	
    alpha <- params[1:nOccCovs]
    gamma <- params[(nOccCovs+1):(nOccCovs+nDetCovs)]	
    eta   <- params[(nOccCovs+nDetCovs+1):(nOccCovs+nDetCovs*2)]	
    
    detection <- Logistic(Xo %*% alpha) * Logistic(Xd %*% gamma) + (1 - Logistic(Xo %*% alpha)) * Logistic(Xd %*% eta)
    return(detection)
}
