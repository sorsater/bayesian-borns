require('msm')
require('mvtnorm')

imgw = 5
imgh = 4
data = read.table('WomenWork.dat.txt', header=TRUE)
Y = data$Work
X = as.matrix(data[,-data$Work])
nObs = nrow(X)

nFeats = ncol(X)
mu_0 = 0
tau = 10
mu_0 = rep(0, nFeats)
omega_0 = tau^2 * diag(nFeats)
beta_prior = as.matrix(rnorm(n=nFeats, mean=mu_0, sd=sqrt(diag(omega_0))))

uGenerator = function(betas) {
  mean = X %*% t(betas)
  u = rep(0, nObs)
  for(i in 1:nObs){
    # Define bounds
    if (Y[1] == 0){
      bounds = c(-Inf, 0)
    }
    else{
      bounds = c(0, Inf)
    }
    u[i] = rtnorm(n=1, mean=mean[i], sd=1, lower=bounds[1], upper=bounds[2])
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
u = uGenerator(t(beta_prior))
draws = 500
result_beta = matrix(0, draws, nFeats)
for(i in 1:draws) {
  print(i)
  beta = betaGenerator(u)
  u = uGenerator(beta)
  result_beta[i,] = beta
}

colnames(result_beta) = colnames(X)
# Posterior mean for beta
res = t(as.matrix(apply(result_beta, 2, mean)))
res


# c
nFeats = ncol(X)
mu = rep(0, nFeats)
tau = 10
sigma2 = tau^2 * diag(nFeats)
beta_prior = rnorm(1, mu, sigma2)


logPostProbit <- function(betas, y, X, mu, sigma2) {
  nPara = length(betas)
  yPred = as.matrix(X) %*% betas
  
  logLike = sum(y*pnorm(yPred, log.p = TRUE) + (1-y)*pnorm(yPred, log.p = TRUE, lower.tail = FALSE))
  
  # evaluating the prior
  logPrior = dmvnorm(betas, matrix(0, nPara, 1), sigma2, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLike + logPrior)
  
}

initValues = rep(0, nFeats)
optimResults = optim(initValues,
                     logPostProbit, y=Y, X=X, mu=mu, sigma2=sigma2,
                     method=c('BFGS'), control=list(fnscale=-1), hessian=TRUE)

postMode = optimResults$par
postCov = -solve(optimResults$hessian)
approxPostStd = sqrt(diag(postCov))

names(postMode) = colnames(X)
names(approxPostStd) = colnames(X)


pdf('plots/betas.pdf')
par(mfrow=c(4,2))
for(i in 1:ncol(X)) {
  betaGrid = seq(min(result_beta[,i]), max(result_beta[,i]), length = 1000)
  approx_density = dnorm(betaGrid, postMode[i], approxPostStd[i])
  gibbs_density = hist(result_beta[,i], breaks=30, plot=FALSE)
  
  plot(gibbs_density,
       #xlim=qnorm(c(0.05, 0.95), postMode[i], approxPostStd[i]),
       ylim=c(0, max(approx_density*sum(gibbs_density$counts), gibbs_density$counts)),
       xlab=colnames(X)[i], col='dodgerblue'
       )
  lines(betaGrid, approx_density*sum(gibbs_density$counts), col='tomato')
}
dev.off()
