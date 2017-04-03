a = 2
b = 2
n = seq(2,2000,10)
p = 0.7
xGrid <- seq(0.001, 0.999, by=0.001)

col1 = '#247ba0'
col2 = '#f25f5c'
hej = function(n){
  normalizedLikelihood = dbeta(xGrid, n*p+1, n*(1-p)+1)
  prior = dbeta(xGrid, a, b)
  posterior = dbeta(xGrid, a+n*p, b+n*(1-p))
  maxDensity <- max(normalizedLikelihood, prior, posterior) # Use to make the y-axis high enough
  #return (c(mean(posterior), sd(posterior)))
  return( c(mean(abs(normalizedLikelihood - posterior)),
            abs(sd(normalizedLikelihood) - sd(posterior))) )
}

data = sapply(n, hej)
plot(n, data[1,], xlab='Samples', ylab='Distance', type='l', lwd=2, col=col1, ylim=c(0, max(data)))
lines(n, data[2,], type='l', lwd=2, col=col2)
legend('topright', c('Mean Error','Standard Deviation'), fill=c(col1, col2), inset=0.02)


n = 10000
p = 0.7
true = pbeta(xGrid, n*p, n*(1-p))
