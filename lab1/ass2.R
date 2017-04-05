#require('geoR')
col1 = '#247ba0'
col2 = '#f25f5c'
imgw = 9
imgh = 6

mu = 3.5
obs = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
n = length(obs)

tausq = sum((log(obs) - mu)^2) / (n)

chiDraw = function(nDraws){
  draws = rchisq(nDraws, n)
  sigmasq = n * tausq / draws
  return (sigmasq)
}

nDraws = 10000
sigmasq = chiDraw(nDraws)

s = seq(0, 1.5, length.out = 1000)
exp_part = exp( ((-n*tausq) / (2*s)) ) / (s^(1+n/2))
constant = ((tausq * n/2)^(n/2)) / (factorial(n/2-1))
CDF = constant * exp_part

pdf('plots/chi-squared.pdf', width=imgw, height=imgh)
  hist(sigmasq, 500, col=col1, border=col1, main='', freq=FALSE, xlim=c(0,1.5))
  # Dependent on geoR
  #lines(s, dinvchisq(s, n, tausq), col=col2, lwd=3)
  # Calculated by hand
  lines(s, CDF, col=col2, lwd=2)
  legend('topright', c('Draws', 'True'), fill=c(col1, col2))
dev.off()

## 2.b Gini
G = function(sigmasq) {
  return (2 * ecdf(sqrt(sigmasq/2)) - 1)
}

#sapply(sigmasq, G)

vals = sqrt(sigmasq/2)
G = 2*pnorm(vals)-1

pdf('plots/g.pdf', width=imgw, height=imgh)
  hist(G, 100, col=col1, border=col1)
dev.off()

