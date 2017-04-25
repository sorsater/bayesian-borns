require('polynom')
require('MASS')
require('geoR')
set.seed(1234567890)
imgw = 7
imgh = 6

data = read.table('TempLinkoping2016.txt', header=TRUE)
n = nrow(data)

# a
#plot(data$time, data$temp, pch=20)

#quad = lm(temp ~ time + I(time^2), data=data)

#coefs = polynomial(c(quad$coefficients))

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

cols = c('dodgerblue', 'firebrick', 'forestgreen')
pdf('plots/posterior.pdf', width=imgw, height=imgh)
  plot(data, pch=20, cex=0.3, xlab='Time of year', ylab='Temperature')
  lines(data$time, predict(polynomial(betas_mean), data$time), col=cols[1], lwd=2)
  lines(data$time, predict(polynomial(ci[1,]), data$time), col=cols[2], lwd=2)
  lines(data$time, predict(polynomial(ci[2,]), data$time), col=cols[3], lwd=2)
  legend('topleft', inset=0.02,
         legend = c('Posterior Mean', '5% interval', '95% interval'),
         fill = cols)
dev.off()
# d
hot_day = data$time[which.max(predict(polynomial(betas_mean), data$time))]
hot_day = as.numeric(betas_mean[2] / (-betas_mean[3]*2))
hot_date = as.Date(hot_day*n, origin="2016-01-01")
hot_date

# e
