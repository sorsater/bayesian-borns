require('manipulate')
set.seed(12345)
col1 = '#247ba0'
col2 = '#f25f5c'
a = 2
b = 2
s = 14
n = 20
xGrid <- seq(0.001, 0.999, by=0.001)

draw = function(nDraws){
  prior = dbeta(xGrid, a, b)
  posterior = rbeta(nDraws, a+s, b+(n-s))
  true = dbeta(xGrid, a+s, b+(n-s))
  
  hist(posterior, 100, prob=TRUE, col=col1, xlim=c(0,1))
  lines(xGrid, true, col=col2, lwd=3)
  #df = data.frame(post=posterior)
  #ggplot(df, aes(x=posterior, y=..density..)) +
  #  geom_histogram(binwidth = 0.005) +
  #  geom_density()

}

manipulate( 
  draw(n), 
  n = slider(1, 100000, step=1, initial=20)
)


##

prob = function(nDraws){
  prior = dbeta(xGrid, a, b)
  posterior = rbeta(nDraws, a+s, b+(n-s))
  true = pbeta(0.4, a+s, b+(n-s))
  
  message((sum(posterior <= 0.4) / length(posterior)) * 100)
  message(round(true*100, 3))
  
}

nDraws = 10000
prob(nDraws)


##

logodds = function(nDraws){
  prior = dbeta(xGrid, a, b)
  posterior = rbeta(nDraws, a+s, b+(n-s))
  phi = log(posterior / (1 - posterior))
  hist(phi, 100, prob=TRUE, col=col1)
  lines(density(phi), lwd=2, col=col2)
}

nDraws = 10000
logodds(nDraws)
