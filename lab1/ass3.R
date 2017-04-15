col1 = '#247ba0'
col2 = '#f25f5c'
imgw = 5
imgh = 4

y = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mu = 2.39

von_misen = function(y, kappa, mu){
  I_zero = besselI(x=kappa, nu=0)
  numerator = exp(kappa * cos(y - mu))
  return (numerator / (2 * pi * I_zero))
}

exponential = function(kappa, lambda=1){
  return (lambda * exp(-lambda*kappa))
}

calculate_posterior = function(y, kappa, mu) {
  prior = exponential(kappa)
  likelihood = prod(von_misen(y, kappa, mu))
  return (prior * likelihood)
}

# To calculate the mode
get_mode = function(vec, kappas) {
  vec = round(vec, 10)
  uniqv = unique(vec)
  mode = uniqv[which.max(tabulate(match(vec, uniqv)))]
  
  mode_idx = which(vec == mode)[1]
  mode_x = kappas[mode_idx]
  
  return (mode_x)
}

kappas = seq(0, 6, 0.001)

posterior = sapply(kappas, calculate_posterior, y=y, mu=mu)
mode = get_mode(posterior, kappas)

pdf("plots/3-posterior-distribution.pdf", width=imgw, height=imgh)
  plot(kappas, posterior, type='l', lwd=2, col=col1,
       xlab=expression(kappa), ylab="Density")
  abline(v=mode, col=col2, lwd=2)
  legend('topright', c('Posterior of k', 'Mode'), fill=c(col1, col2), inset=0.02)
dev.off()
