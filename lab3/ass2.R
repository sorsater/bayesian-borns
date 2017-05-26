require('msm')
require('mvtnorm')

data = read.table('WomenWork.dat', header=TRUE)
Y = data$Work
X = as.matrix(data[,-data$Work])
nObs = nrow(X)
nFeats = ncol(X)

# Calculate the accuracy of the beta values
performance = function(betas){
  y_hat = as.vector(X %*% betas)
  y_hat = sign(y_hat)
  y_hat[y_hat==-1] = 0
  
  print(table(y_hat, Y))
  print(mean(y_hat != Y))
}

uGenerator = function(betas) {
  curMean = X %*% t(betas)
  u = rep(0, nObs)
  for(i in 1:nObs){
    # Define bounds = (lower, upper)
    if(Y[i] == 0){
      bounds = c(-Inf, 0)
    }else{
      bounds = c(0, Inf)
    }
    u[i] = rtnorm(n=1, mean=curMean[i], sd=1, lower=bounds[1], upper=bounds[2])
  }
  return(u)
}

# From the given u, draw new beta-values
betaGenerator = function(u, sigma2=1) {
  # From lecture 5, slide 8
  XtX = t(X) %*% X
  beta_hat = solve(XtX) %*% t(X) %*% as.matrix(u)
  mu_n = solve(XtX + omega_0) %*% (XtX %*% beta_hat + omega_0 %*% mu_0)
  omega_n = XtX + omega_0
  
  beta = rmvnorm(n=1, mean=mu_n, sigma=(sigma2*solve(omega_n)))
  return (beta)
}

# Initial values
mu_0 = 0
tau = 10
mu_0 = rep(0, nFeats)
covar_0 = tau^2 * diag(nFeats)
omega_0 = solve(covar_0)

beta_prior = as.matrix(rnorm(n=nFeats, mean=mu_0, sd=sqrt(diag(omega_0))))
u = uGenerator(t(beta_prior))

# Necessary with 1000 iterations, draws more to get smoother histograms
draws = 2500
# One column for each beta parameter
result_beta = matrix(0, draws, nFeats)
for(i in 1:draws) {
  print(i)
  beta = betaGenerator(u)
  u = uGenerator(beta)
  result_beta[i,] = beta
}

performance(colMeans(result_beta))

# c

logPostProbit = function(betas, y, X) {
  yPred = as.matrix(X) %*% betas
  
  logLike = sum(y*pnorm(yPred, log.p = TRUE) + (1-y)*pnorm(yPred, log.p = TRUE, lower.tail = FALSE))
  logPrior = dmvnorm(betas, mu_0, covar_0, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLike + logPrior)
  
}

initValues = rep(0, nFeats)
optimResults = optim(initValues,
                     logPostProbit, y=Y, X=X,
                     method=c('BFGS'), control=list(fnscale=-1), hessian=TRUE)

postMode = optimResults$par
postCov = -solve(optimResults$hessian)
approxPostStd = sqrt(diag(postCov))

performance(postMode)

pdf('plots/betas.pdf')
  par(mfrow=c(4,2))
  for(i in 1:nFeats) {
    mean = postMode[i]
    sd = approxPostStd[i]
  
    # Remove outliers. (2 < x < 98) %
    draws = result_beta[,i]
    threshold = 0.02
    bounds = quantile(draws, probs=c(threshold, 1 - threshold))
    draws = draws[draws > bounds[1]]
    trimmed = draws[draws < bounds[2]]
    
    # Toggle to include/exclude outliers
    values = result_beta[,i]
    values = trimmed
    
    xRange = c(values, mean + sd * 4, mean - sd * 4)
    betaGrid = seq(min(xRange), max(xRange), length=1000)
    
    gibbs_density = hist(values, breaks=50, plot=FALSE)
    approx_density = dnorm(betaGrid, mean=mean, sd=sd)
    
    yMax = max(c(approx_density, gibbs_density$density))
    
    plot(gibbs_density, freq=FALSE, col='dodgerblue', border='dodgerblue',
         xlim=c(min(betaGrid), max(betaGrid)), 
         ylim=c(0, yMax),
         xlab=substitute(beta[idx], list(idx=i)),
         main=colnames(X)[i])
    lines(betaGrid, approx_density, col='tomato', lwd=2)
  }
dev.off()

