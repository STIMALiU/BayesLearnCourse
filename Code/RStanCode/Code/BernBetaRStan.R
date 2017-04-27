# BRugs code for analyzing a simple Bernoulli model with a Beta(a,b) prior
# Author: Mattias Villani, Linköping University
# Ported to RStan by Måns Magnusson
# Date:   2013-10-24

rm(list=ls())

#install.packages("rstan")
library(rstan)

# Data 
x = c(1,1,0,0,1,1,1,1,0,1)
n = length(x)
a = 1
b = 1

BernBetaData <- list(n = length(x), x=x, a=1, b=1)


# Model
BernBetaStanModel <- '
data {
  int<lower=0> n;
  int<lower=0,upper=1> x[n];
  real<lower=0> a;
  real<lower=0> b;
}

parameters {
  real<lower=0,upper=1> theta;
} 

model {
  theta ~ beta(a,b);
  for (i in 1:n)
    x[i] ~ bernoulli(theta);
}
'

# Do the fitting of the model
fit1<-stan(model_code=BernBetaStanModel,
           data=BernBetaData,
           warmup=1000,
           iter=2000,
           chains=4)

print(fit1,digits_summary=3)
plot(fit1)

# Plot some results
res <- extract(fit1) 
thetaSeq <- seq(0,1,by=0.01)
par(mfrow = c(1,1))
hist(res$theta, 40, freq = FALSE, main = 'Posterior of theta - all chains', xlab ='theta') # histogram of draws of a from the first Markov Chain
lines(thetaSeq, dbeta(thetaSeq, shape1 = sum(x) + a, shape2 = n - sum(x) + b), col = "red")
legend("topleft", inset=.05, legend = c('MCMC approximation','True density'), lty =c(1,1),col=c('black','red'))

