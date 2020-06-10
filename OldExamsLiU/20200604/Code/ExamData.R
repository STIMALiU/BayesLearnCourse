# Bayesian Learning Exam 2020-06-04
# Run this file once during the exam to get all the required data and functions for the exam in working memory 
# Author: Per Siden

###############################
########## Problem 1 ########## 
############################### 

###############################
########## Problem 2 ########## 
############################### 

###############################
########## Problem 3 ########## 
############################### 

# Reading the data from file
load(file = 'titanic.RData')

if("mvtnorm" %in% rownames(installed.packages()) == FALSE) {install.packages("mvtnorm")}
if("msm" %in% rownames(installed.packages()) == FALSE) {install.packages("msm")}
library(mvtnorm) # For mulitvariate normal
library(msm) # For truncated normal

BayesProbReg <- function(y, X, mu_0, tau, nIter){
  # Gibbs sampling in probit regression using data augmentation:
  #
  # beta | tau ~ N(mu_0, tau^2*I)
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   tau - prior standard deviation for beta
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   betaSample    - Posterior samples of beta.          nIter-by-nCovs matrix
  
  # Prior
  nPara <- dim(X)[2] # this line was missing in the exam
  priorCov <- tau^2*diag(nPara)
  priorPrec <- solve(priorCov)
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  n1 = sum(y)
  n0 = n - n1
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  
  # The actual sampling
  betaSample = matrix(NA, nIter, nCovs)
  u <- matrix(NA, n, 1)
  beta <- solve(XX,crossprod(X,y)) # OLS estimate as initial value
  for (i in 1:nIter){
    
    xBeta <- X%*%beta
    
    # Draw u | beta
    u[y == 0] <- rtnorm(n = n0, mean = xBeta[y==0], sd = 1, lower = -Inf, upper = 0)
    u[y == 1] <- rtnorm(n = n1, mean = xBeta[y==1], sd = 1, lower = 0, upper = Inf)
    
    # Draw beta | u
    betaHat <- solve(XX,t(X)%*%u)
    postPrec <- XX + priorPrec
    postCov <- solve(postPrec)
    betaMean <- solve(postPrec,XX%*%betaHat + priorPrec%*%mu_0)
    beta <- t(rmvnorm(n = 1, mean = betaMean, sigma = postCov))
    betaSample[i,] <- t(beta)
    
  }
  
  return(betaSample=betaSample)
}


###############################
########## Problem 4 ########## 
###############################