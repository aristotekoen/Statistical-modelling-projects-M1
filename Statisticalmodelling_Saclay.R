rm(list=objects())
graphics.off()
B=200
sigma=0.2
gamma=seq(0.1,0.9,0.1)
N=c(20,50,100)
foo= function(x) x*sigma*sigma/(1-x)
curve(foo(x),ylim=c(0,2))


RiskMLE = function(gamma, n){
  d = floor(gamma*n)
  theta = rep(0.1, d)
  X = matrix(rnorm(d*n), nrow = n, ncol = d)
  eps = rnorm(n, 0, sigma)
  Y = X%*%theta + eps
  theta_hat = solve(t(X)%*%X) %*% t(X) %*% Y 
  r = t(theta_hat-theta)%*%(theta_hat-theta)
  return(r)
}

col = 1
for (i in N){
  Risk = 1:length(gamma)
  k = 1
  for (j in gamma){
    R = replicate(B, RiskMLE(j, i))
    Risk[k] = mean(R)
    k = k+ 1
  }
  col = col + 1
  points(gamma, Risk, type = 'p', xlim = c(0,1), pch=16, col = col)
}

legend('topleft',c('n=20','n=50','n=100'),col = 2:4, pch = 16)
legend('topright',"risque asymptotique" ,col = 'black', lty = 1)

#Q13

n = 100
d = 30
N = 20
lambda_ = 0:N

X = matrix(rnorm(d*n),nrow = n,ncol = d)
theta = rep(0.1, d)
eps = rnorm(n, 0, sigma)
Y = X%*%theta + eps

bia = 1:B; bia.est = 1:N+1
var = 1:B; var.est = 1:N+1
ris = 1:B; ris.est = 1:N+1

for (lambda in lambda_){
  for (k in 1:B){
    theta.chap = 1/n*solve(1/n*t(X)%*%X + lambda*diag(d))%*%t(X)%*%Y
    theta.barr = mean(theta.chap)
    bia[k] = t(theta.barr - theta     )%*%(theta.barr - theta     )
    var[k] = t(theta.chap - theta.barr)%*%(theta.chap - theta.barr)
    ris[k] = t(theta.chap - theta     )%*%(theta.chap - theta     )
  }
  bia.est[lambda+1] = mean(bia)
  var.est[lambda+1] = mean(var)
  ris.est[lambda+1] = mean(ris)
}

par(mfrow = c(1,3))
plot(lambda_, bia.est,
     ylab = "biais estimÈ", pch="+", col="red")
plot(lambda_, var.est,
     ylab = "variance estimÈe", pch="+", col="blue")
plot(lambda_, ris.est,
     ylab = "risque estimÈ", pch="+", col="black")
title("Estimation des propriÈtÈes du modËle")



#Q14


B=200
sigma=0.2
lamda=4
gamma=seq(0.1,0.9,0.1)
N=c(20,50,100)
foo= function(x) x*sigma*sigma/(1-x)
curve(foo(x),ylim=c(0,2))


RiskLamda = function(gamma, n){
  d = floor(gamma*n)
  theta = rep(0.1, d)
  X = matrix(rnorm(d*n), nrow = n, ncol = d)
  eps = rnorm(n, 0, sigma)
  Y = X%*%theta + eps
  theta_hat_lamda = 1/n * solve(1/n * t(X)%*%X + lamda * diag(d)) %*% t(X)%*%Y
  r = t(theta_hat_lamda-theta)%*%(theta_hat_lamda-theta)
  return(r)
}

col = 1
for (i in N){
  Risk = 1:length(gamma)
  k = 1
  for (j in gamma){
    R = replicate(B, RiskLamda(j, i))
    Risk[k] = mean(R)
    k = k+ 1
  }
  col = col + 1
  points(gamma, Risk, type = 'p', xlim = c(0,1), pch=16, col = col)
}

legend('topleft',c('n=20','n=50','n=100'),col = 2:4, pch = 16)
legend('topright',"risque asymptotique" ,col = 'black', lty = 1)



