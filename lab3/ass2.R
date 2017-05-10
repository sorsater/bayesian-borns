require('msm')
imgw = 5
imgh = 4
data = read.table('WomenWork.dat.txt', header=TRUE)
Y = data$Work
X = as.matrix(data[,-data$Work])


nFeats = ncol(X)
mu = 0
tau = 10
beta_prior = rnorm(n=nFeats, mean=mu, sd=tau)
plot(density(beta_prior))


betaGenerator = function(u, y) {
  beta_hat = solve(t(u) %*% u) %*% t(u) %*% y
  return (beta_hat)
}

uGenerator = function(x, beta, sigma2=1) {
  mean = t(x)%*%beta
  sd = sqrt(sigma2)
  print(mean)
  draw = rnorm(n=1, mean=mean, sd=sd)
  
  return(draw)
}

#u = 
#u = c()
for(i in 1:10) {
  u = c(u, uGenerator(X, beta))
  beta = c(beta, betaGenerator(u, Y))
}

# stopp























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
  