require('MASS')
require('geoR')
require('mvtnorm')

ebay = read.table('eBayNumberOfBidderData.dat', header = TRUE)
X = as.matrix(ebay[,-1])
Y = ebay[,1]

nFeats = ncol(X)
nObs = nrow(X)

# a
glmModel = glm(nBids ~ 0 + ., data=ebay, family=poisson)
coefs = glmModel$coefficients
coefs

# b
sigma2_0 = 100 * solve(t(X) %*% X)
mu_0 = rep(0, nFeats)

beta_prior = mvrnorm(n=1, mu=mu_0, Sigma=sigma2_0)


logPostPoi = function(betas, y, X) {
  yPred = as.matrix(X) %*% betas
  
  lambda = exp(yPred)
  #logLike = sum(log(exp(-lambda) * lambda^y / factorial(y)))
  logLike = sum(log(exp(-lambda) * lambda^y / factorial(y)))
  logPrior = dmvnorm(betas, mu_0, sigma2_0, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLike + logPrior)
  
}

initValues = rep(0, nFeats)
optimResults = optim(initValues,
                     logPostPoi, y=Y, X=X,
                     method=c('BFGS'), control=list(fnscale=-1), hessian=TRUE)

postMode = optimResults$par
postCov = -solve(optimResults$hessian)
approxPostStd = sqrt(diag(postCov)/nFeats)
postMode

# c


metroDraw = function(logPostFunc, thetas, sigma, c_tilde, ...) {
  draws = 2000
  result_theta = matrix(0, draws, nFeats)
  mean_conv = matrix(0, draws, nFeats)
  thetas_c = thetas
  #prob = c()
  for(i in 1:draws){
    thetas_p = (mvrnorm(1, thetas, c_tilde*sigma))
    alpha = min(1, exp(logPostFunc(thetas_p, ...) - logPostFunc(thetas_c, ...)))
    
    update = sample(c(FALSE, TRUE), 1, prob=c(1-alpha, alpha))
    #prob = c(prob, alpha)
    if(update)
      thetas_c = thetas_p
    
    result_theta[i,] = thetas_c
    
    if (i > 1)
      mean_conv[i,] = colMeans(result_theta[1:i,])
    
      
  }
  #print(mean(prob))
  return (mean_conv)
}


c_tilde = 0.5
thetas = metroDraw(logPostPoi, thetas=postMode, sigma=postCov, c_tilde=c_tilde, y=Y, X=X)
#plot(thetas[,1], type='l')
par(mfrow=c(3,3))
for(i in 1:nFeats) {
  plot(thetas[,i], type='l', col=rainbow(nFeats)[i], xlab='Samples', ylab=colnames(X)[i])
}

# d
optTheta = thetas[nrow(thetas),]
femkrona_1972 = c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)
names(femkrona_1972) = colnames(X)
lambda = exp(femkrona_1972%*%optTheta)
xGrid = 0:5
plot(xGrid, ppois(xGrid, lambda), cex=1, pch=20, xlab='k', ylab='P(x<=k)', col='dodgerblue')
text(0.2, 0.4, paste(100*round(ppois(0, lambda), 2), '%'), col='tomato')


