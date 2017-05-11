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
  beta = betaGenerator(u)cd
  u = uGenerator(beta)
  result_beta[i,] = beta
}

# Posterior mean for beta
res = t(as.matrix(apply(result_beta, 2, mean)))
colnames(res) = colnames(X)
res