col1 = '#247ba0'
col2 = '#f25f5c'
imgw = 9
imgh = 6

obs = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
mu = 2.39

von_misen = function(y, kappa, mu){
  I_zero = besselI(x=kappa, nu=0)
  numerator = exp(kappa * cos(y - mu))
  return (numerator / (2 * pi * I_zero))
}

exponential = function(kappa, lambda=1){
  return (lambda * exp(-lambda*kappa))
}

# 
calculate_mode = function(vec) {
  uniqv = unique(vec)
  uniqv[which.max(tabulate(match(vec, uniqv)))]
}

kappas = seq(0, 6, 0.001)
res = c()
for (kappa in kappas){
  likelihood = prod(von_misen(obs, kappa, mu))
  prior = exponential(kappa)
  
  #res = c(res, likelihood)
  res = c(res, likelihood * prior)
}

# Need to tune the rounding parameter, all values are unique and need to trim off some decimals
rounded_res = round(res,10)
mode_idx = which(rounded_res == calculate_mode(rounded_res))[1]
mode = kappas[mode_idx]

pdf("plots/posterior-distribution.pdf", width=imgw, height=imgh)
  plot(kappas, res, type='l', lwd=2, col=col1, xlab="Kappa", ylab="Density")
  abline(v=mode, col=col2, lwd=2)
  legend('topright', c('Posterior of k', 'Mode'), fill=c(col1, col2), inset=0.02)
dev.off()

