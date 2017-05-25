set.seed(12345)
col1 = '#247ba0'
col2 = '#f25f5c'
col3 = '#333333'

imgw = 5
imgh = 4

a = 2
b = 2
s = 14
n = 20

xGrid <- seq(0.001, 0.999, 0.001)

## 1a
draw = function(nDraws){
  sample = rbeta(nDraws, a+s, b+(n-s))
  mean_sample = mean(sample)
  sd_sample = sd(sample)
  
  return (c(mean_sample, sd_sample))
}

alpha = a + s
beta = b + (n-s)

intervals = seq(20,50000,100)
data = sapply(intervals, draw)

true = dbeta(xGrid, a+s, b+(n-s))
mean_true = alpha / (alpha + beta)
sd_true = sqrt( (alpha*beta) / ((alpha+beta)^2*(alpha+beta+1)) )

# Mean
pdf("plots/1-converge_mean.pdf", width=imgw, height=imgh)
  plot(intervals, abs(data[1,]-mean_true), type='l', col=col3,
    xlab = 'nDraws', ylab='Absolute error')
dev.off()

# SD
pdf("plots/1-converge_sd.pdf", width=imgw, height=imgh)
  plot(intervals, abs(data[2,]-sd_true), type='l', col=col3, 
    xlab = 'nDraws', ylab='Absolute error')
dev.off()

## 1b
nDraws = 10000
posterior = rbeta(nDraws, a+s, b+(n-s))
true = pbeta(0.4, a+s, b+(n-s))

# Simulated value
message((sum(posterior <= 0.4) / length(posterior)) * 100)
# Theoretical value
message(round(true*100, 3))

## 1c
nDraws = 10000
posterior = rbeta(nDraws, a+s, b+(n-s))
phi = log(posterior / (1 - posterior))

pdf("plots/1-phi.pdf", width=imgw, height=imgh)
  hist(phi, 100, prob=TRUE, col=col1, border=col1, main='', xlab=expression(phi))
  lines(density(phi), lwd=2, col=col2)
dev.off()
