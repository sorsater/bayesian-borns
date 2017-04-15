col1 = '#247ba0'
col2 = '#f25f5c'
col3 = '#333333'
imgw = 5
imgh = 4

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

pdf('plots/2-chi-squared.pdf', width=imgw, height=imgh)
  hist(sigmasq, 500, col=col1, border=col1, main='', freq=FALSE, xlim=c(0,1.5))
  lines(s, CDF, col=col2, lwd=2)
  legend('topright', c('Draws', 'True'), fill=c(col1, col2))
dev.off()

## 2.b Gini
vals = sqrt(sigmasq/2)
G = 2 * pnorm(vals) - 1

pdf('plots/2-g.pdf', width=imgw, height=imgh)
  hist(G, 100, col=col1, border=col1, main='')
dev.off()


## 2.c Credible Intervals

# Equal tail interval
x_eti = quantile(G, probs=c(0.025, 0.975))

# Highest Posterior Density
dg = density(G)
ordered_x = dg$x[order(-dg$y)]
ordered_y = dg$y[order(-dg$y)]
  
cur = 0
for(i in 1:length(dg$y)){
  cur = cur + ordered_y[i]
  if(cur / sum(dg$y) >= 0.95){
    break
  }
}
x_hpd = c(min(ordered_x[1:i]), max(ordered_x[1:i]))

pdf('plots/2-credible-intervals.pdf', width=imgw, height=imgh)
  plot(dg, col=col3, lwd=2, main="")
  lines(x_eti, rep(0.12, 2), col=col1, lwd=4)
  lines(x_hpd, rep(0.02, 2), col=col2, lwd=4)
  legend("topright", 
         legend = c("Equal Tail Interval","Highest Posterior Density"),
         fill = c(col1, col2),
         inset = 0.02)
dev.off()


