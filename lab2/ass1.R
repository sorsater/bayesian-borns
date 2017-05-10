require('polynom')
require('MASS')
require('geoR')
set.seed(1234567890)
imgw = 5
imgh = 4

data = read.table('TempLinkoping2016.txt', header=TRUE)
n = nrow(data)

# a
quad = lm(temp ~ time + I(time^2), data=data)
coefs = polynomial(c(quad$coefficients))

#plot(data$time, data$temp, pch=20)
#lines(data$time, predict(coefs, data$time))

# Hyper parameters
mu_0 = c(-10, 100, -100)
omega_0 = diag(3)
v_0 = n - 3
sigma2_0 = 5

# b
cols = rainbow(5)

pdf('plots/hyper.pdf', width=imgw, height=imgh)
  plot(data, pch=20, cex=0.3, type='n', xlab='Time of year', ylab='Temperature')
  for (i in 1:5){
    sigma2 = rinvchisq(1, df=v_0, scale=sigma2_0)
    betas = mvrnorm(n = 1, mu = mu_0, Sigma=sigma2 * solve(omega_0))
    coefs = polynomial(betas)
    lines(data$time, predict(coefs, data$time), col=cols[i])
  }
dev.off()

# c
y = data$temp
X = as.matrix(data.frame(beta0=rep(1, n), beta1=data$time, beta2=data$time^2))

beta_hat = solve(t(X) %*% X) %*% t(X) %*% y
mu_n = solve(t(X) %*% X + omega_0) %*% (t(X) %*% X %*% beta_hat + omega_0 %*% mu_0)
omega_n = t(X) %*% X + omega_0
v_n = v_0 + n
v_n.sigma2_n = v_0 * sigma2_0 +
                (t(y) %*% y +
                t(mu_0) %*% omega_0 %*% mu_0 -
                t(mu_n) %*% omega_n %*% mu_n)

sigma2 = rinvchisq(1000, df=v_n, scale=v_n.sigma2_n/v_n) 
betas = sapply(sigma2, function(sigma2){mvrnorm(n=1, mu=mu_n, Sigma=sigma2*solve(omega_n))})
betas_mean = apply(betas, 1, mean)
print(betas_mean)

ci = apply(betas, 1, quantile, probs=c(0.05, 0.95))

cols = c('forestgreen', 'dodgerblue', 'tomato')
pdf('plots/posterior.pdf', width=imgw, height=imgh)
  plot(data, pch=20, cex=0.3, xlab='Time of year', ylab='Temperature', xaxt='n')
  lines(data$time, predict(polynomial(ci[2,]), data$time), col=cols[1], lwd=2)
  lines(data$time, predict(polynomial(betas_mean), data$time), col=cols[2], lwd=2)
  lines(data$time, predict(polynomial(ci[1,]), data$time), col=cols[3], lwd=2)
  legend('topleft', inset=0.02,
         legend = c('95% interval', 'Posterior Mean','5% interval'),
         fill = cols, cex=0.65)
  axis(1, at=c(seq(0,1,length=12)), 
       labels=c('Jan', 'Feb','Mar', 'Apr','May', 'Jun','Jul', 'Aug', 'Sep', 'Oct','Nov', 'Dec'), cex=0.4)
dev.off()
# d
hot_day = data$time[which.max(predict(polynomial(betas_mean), data$time))]
hot_day
hot_day = as.numeric(betas_mean[2] / (-betas_mean[3]*2))
hot_day
hot_date = as.Date(hot_day*n, origin="2016-01-01")
hot_date

xGrid = seq(0,1,0.001)
a = dnorm(xGrid, mean=hot_day, sd=0.1)
tjo = max(a)
plot(as.Date(xGrid*n, origin="2016-01-01"),a, type='l')
abline(h=tjo/2)
