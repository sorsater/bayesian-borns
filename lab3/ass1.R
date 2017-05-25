require('geoR')

imgw = 7
imgh = 4
rainfall = read.table('rainfall.dat')[,1]
set.seed(12345)
n = length(rainfall)

# a
mu_0 = mean(rainfall)
sigma2_0 = 1
v_0 = 1

nDraws = 2500

sigma2 = rinvchisq(n=1, v_0, sigma2_0)
# Lecture 7 slide 17
result = matrix(0, nDraws, 4)
for(i in 1:nDraws) {
  mu = rnorm(n=1, mean=mu_0, sd=sqrt(sigma2 / n))
  sigma2 = rinvchisq(n=1, df=(v_0 + n), scale=(v_0*sigma2_0 +  sum((rainfall - mu)^2)) / (n + v_0))
  result[i,] = c(mu, sigma2, mean(result[1:i-1,1]), mean(result[1:i-1,2]))
}

#pdf('plots/mu.pdf', width=imgw, height=imgh)
#plot(result[,3], type='l', xlab='Samples', ylab='Mean', col='dodgerblue')
#dev.off()
#pdf('plots/sigma.pdf', width=imgw, height=imgh)
#plot(result[,4], type='l', xlab='Samples', ylab='Variance', col='tomato')
#dev.off()

# Plot mean and variance in same plot
d = data.frame(x=seq(1, nDraws), mu=result[,3], sigma=result[,4])

pdf('plots/mean_variance.pdf', width=imgw, heigh=imgh)
par(mar = c(5,5,2,5))
with(d, plot(x, mu, type='l', col='dodgerblue', ylab='Mean', lwd=2, xlab='Samples'))
par(new = T)
with(d, plot(x, sigma, type='l', col='tomato', axes=F, xlab=NA, ylab=NA, lwd=2))
axis(side = 4)
mtext(side=4, line=3, 'Variance')
legend('right', legend=c('Mean', 'Variance'), fill=c('dodgerblue', 'tomato'), inset=0.02)
dev.off()

# b

# Using template for this part

# Estimating a simple mixture of normals
# Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com

x = as.matrix(read.table('rainfall.dat')[,1])

# Model options
nComp = 2    # Number of mixture components

# Prior options
alpha = 10*rep(1,nComp) # Dirichlet(alpha)
muPrior = c(x[which.max(density(x)$y)], mean(x)) # Prior mean of theta
tau2Prior = rep(10,nComp) # Prior std theta
sigma2_0 = rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 = rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter = 100 # Number of Gibbs sampling draws

# Plotting options
plotFit = TRUE
lineColors = c("dodgerblue", "lawngreen", "magenta", 'yellow')
sleepTime = 0.02 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 = function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet = function(param){
  nCat = length(param)
  thetaDraws = matrix(NA,nCat,1)
  for (j in 1:nCat){
    thetaDraws[j] = rgamma(1,param[j],1)
  }
  # Diving every column of ThetaDraws by the sum of the elements in that column.
  thetaDraws = thetaDraws/sum(thetaDraws) 
  return(thetaDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc = function(S){
  n = dim(S)[1]
  alloc = rep(0,n)
  for (i in 1:n){
    alloc[i] = which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs = length(x)
# nObs-by-nComp matrix with component allocations.
S = t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) 
theta = quantile(x, probs = seq(0,1,length = nComp))
sigma2 = rep(var(x),nComp)
probObsInComp = rep(NA, nComp)

# Setting up the plot
xGrid = seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 500)
xGridMin = min(xGrid)
xGridMax = max(xGrid)
mixDensMean = rep(0,length(xGrid))
effIterCount = 0
ylim = c(0,1.33*max(hist(x, breaks=50)$density))

result_mu = matrix(0, nIter, 4)
result_sigma = matrix(0, nIter, 4)
for (k in 1:nIter){
  message(paste('Iteration number:',k))
  # Just a function that converts between different representations of the group allocations
  alloc = S2alloc(S) 
  nAlloc = colSums(S)
  print(nAlloc)
  # Update components probabilities
  w = rDirichlet(alpha + nAlloc)
  
  # Update theta's
  for (j in 1:nComp){
    precPrior = 1/tau2Prior[j]
    precData = nAlloc[j]/sigma2[j]
    precPost = precPrior + precData
    wPrior = precPrior/precPost
    muPost = wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post = 1/precPost
    theta[j] = rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  result_mu[k,] = c(theta[1], theta[2], mean(result_mu[1:k-1,1]), mean(result_mu[1:k-1,2]))
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] = rScaledInvChi2(1, df = nu0[j] + nAlloc[j],
              scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - theta[j])^2))/(nu0[j] + nAlloc[j]))
  }
  result_sigma[k,] = c(sigma2[1], sigma2[2], mean(result_sigma[1:k-1,1]), mean(result_sigma[1:k-1,2]))
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] = w[j]*dnorm(x[i], mean = theta[j], sd = sqrt(sigma2[j]))
    }
    S[i,] = t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount = effIterCount + 1
    hist(x, breaks = 50, freq = FALSE,
         xlim = c(xGridMin,250),
         ylim = ylim, col="gainsboro",
         main = paste("Iteration number",k))
    mixDens = rep(0,length(xGrid))
    components = c()
    for (j in 1:nComp){
      compDens = dnorm(xGrid,theta[j],sd = sqrt(sigma2[j]))
      mixDens = mixDens + w[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] = paste("Component ",j)
    }
    mixDensMean = ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'tomato')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'), 
           col = c("black",lineColors[1:nComp], 'tomato'), lwd = 2)
    #Sys.sleep(sleepTime)
  }
  
}
pdf('plots/mixture.pdf', width=7, height=imgh)
  hist(x, breaks = 100, freq = FALSE, 
     xlim = c(xGridMin,250), main = "",
     col='gainsboro',
     xlab='Precipitation')
  lines(xGrid, mixDensMean, type = "l", lwd = 2, col = "tomato")
  lines(xGrid, 
      dnorm(xGrid, mean = result[nrow(result),3], sd = sqrt(result[nrow(result),4])), 
      type = "l", lwd = 2, col = "dodgerblue")
  legend("topright", box.lty = 1, 
       legend = c("Data histogram","Normal density (a)", "Mixture density (b)"), 
       col=c("gainsboro","dodgerblue","tomato"), lwd = 2)
dev.off()

ylim_mu = c(min(result_mu[-1,3:4]), max(result_mu[-1,3:4]))
ylim_sigma = c(min(result_sigma[-1,3:4]), max(result_sigma[-1,3:4]))

# Plot the result
pdf('plots/muMixed.pdf', width=5, height=imgh)
  plot(result_mu[,3], type='l', xlab='Samples', ylab='Mean', col='dodgerblue',
     ylim=ylim_mu, lwd=2)
  lines(result_mu[,4], col='tomato', lwd=2)
  legend("right", box.lty = 1, inset=0.02,
       legend=c("Component 1","Component 2"), 
       col=c("dodgerblue", "tomato"), lwd = 2)
dev.off()


# Plot the result but with traceplot
mcmc_mu = mcmc(result_mu[,1:2])
ylims = c(min(mcmc_mu), max(mcmc_mu))
pdf('plots/muMixed-trace.pdf', width=5, height=imgh)
traceplot(mcmc_mu[,1], type='l', xlab='Samples', ylab='Mean', col='dodgerblue',
     ylim=ylims,#ylim_mu,
     lwd=2)
lines(mcmc_mu[,2], col='tomato', lwd=2)
legend("right", box.lty = 1, inset=0.02,
       legend=c("Component 1","Component 2"), 
       col=c("dodgerblue", "tomato"), lwd = 2)
dev.off()



pdf('plots/sigmaMixed.pdf', width=5, height=imgh)
  plot(result_sigma[,3], type='l', xlab='Samples', ylab='Variance', col='dodgerblue',
     ylim=ylim_sigma, lwd=2)
  lines(result_sigma[,4], col='tomato', lwd=2)
  legend("right", box.lty = 1, inset=0.02,
    legend=c("Component 1","Component 2"), 
    col=c("dodgerblue", "tomato"), lwd = 2)

dev.off()

