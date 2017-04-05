require('geoR')
col1 = '#247ba0'
col2 = '#f25f5c'

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

hist(sigmasq, 500, col=col1, border=col1, main='', freq=FALSE, xlim=c(0,1.5))
lines(s, dinvchisq(s, n, tausq), col=col2, lwd=3)

#PDF = exp( ((-n*tausq) / (2*s)) ) / (s^(6))

## 2.b Gini
G = function(sigmasq) {
  return (2 * ecdf(sqrt(sigmasq/2)) - 1)
}

sapply(sigmasq, G)


vals = sqrt(sigmasq/2)
G = 2*pnorm(vals)-1
hist(G, 100)

getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mode = getmode(G)
mode
