require('mvtnorm')
require('LaplacesDemon')
imgw = 5
imgh = 4
data = read.table('WomenWork.dat.txt', header=TRUE)
response = data$Work
features = as.matrix(data[,-data$Work])

# a
fit = glm(Work ~ 0 + ., data = data, family = 'binomial')

# b
nFeats = ncol(features)
mu = rep(0, nFeats)
tau = 10
sigma2 = tau^2 * diag(nFeats)
beta_prior = rnorm(mu, sigma2)

logPostLogistic = function(betas, y, X, mu, sigma2) {
  nPara = length(betas)
  yPred = as.matrix(X) %*% betas
  
  logLike = sum(yPred*y - log(1 + exp(yPred)))
  if (abs(logLike) == Inf) {
    logLike = -20000
  }
  logPrior = dmvnorm(betas, rep(0, nPara), sigma2, log=TRUE)
  logPost = logPrior + logLike
  return (logPost)
}

initValues = rep(0, nFeats)
optimResults = optim(initValues,
                     logPostLogistic, y=response, X=features, mu=mu, sigma2=sigma2,
                     method=c('BFGS'), control=list(fnscale=-1), hessian=TRUE)

postMode = optimResults$par
postCov = -solve(optimResults$hessian)
approxPostStd = sqrt(diag(postCov))

names(postMode) = colnames(features)
names(approxPostStd) = colnames(features)

childMean = postMode['NSmallChild']
childStd = approxPostStd['NSmallChild']

betaGrid = seq(-abs(childMean - 4*childStd), abs(childMean + 4* childStd), length = 1000)

pdf('plots/smallchildren.pdf', width=imgw, height=imgh)
  plot(betaGrid, dnorm(betaGrid, childMean, childStd), 
       type = "l", lwd = 2, ylab = '', xlab = expression(beta), col='steelblue')
  lines(qnorm(c(0.025, 0.975), childMean, childStd), c(0,0), lwd=2, col='tomato')
  legend('topright', legend=c('95% equal tail interval'), fill=c('tomato'), inset=0.02, cex=0.6)
dev.off()
# c

inverseLogit = function(betas, features){
  return (exp(features %*% betas) / (1 + exp(features %*% betas)))
}

jane_doe = c(1, 10, 8, 10, 1, 40, 1, 1)
betas = rmvnorm(1000, mean=postMode, sigma=postCov)

p = apply(betas, 1, inverseLogit, features=jane_doe)
y = sapply(p, rbern, n=1)

h = hist(y, breaks=2, plot=FALSE)
h$counts = h$counts / sum(h$counts)
h$counts = h$counts * 100
names(h$counts) = h$counts

pdf('plots/working.pdf', width=imgw, height=imgh)
  barplot(h$counts, ylab='Confidence', ylim=c(0,100), col=c('tomato', 'dodgerblue'), 
          main='', legend=c('Not working', 'Working'))
dev.off()

