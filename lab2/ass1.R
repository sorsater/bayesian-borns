require('polynom')
require('MASS')
require('geoR')

data = read.table('TempLinkoping2016.txt', header=TRUE)
n = nrow(data)

# a
#plot(data$time, data$temp, pch=20)

#quad = lm(temp ~ time + I(time^2), data=data)

#coefs = polynomial(c(quad$coefficients))

#plot(data$time, data$temp, pch=20)
#lines(data$time, predict(coefs, data$time))

# Hyper parameters
mu_0 = c(-10, 40, -40)
omega_0 = diag(3)
v_0 = n - 3
sigma2_0 = 1

# b

for (i in 1:5){
  betas = mvrnorm(n = 1, mu = mu_0, Sigma=sigma2_0 * solve(omega_0))
  sigma2 = rchisq(n = v_0, )
  
}


