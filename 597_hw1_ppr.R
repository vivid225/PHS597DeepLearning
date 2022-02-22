#  Implement project pursuit algorithm and compare it with ppr package in R. 


## Simulate data -----

n=20

sigma<-diag(c(1,1))
X<-mvrnorm(n,c(0,3),sigma)
x1 <- X[,1]
x2 <- X[,2]
y <- x1+log(x1^2) + exp(x1)


## ppr results -----
pprfit <- ppr(X,y,nterms = 1, max.terms = 1)
summary(pprfit)

## Build functions -----

### initiate beta
b<-c(1/2,1/2)

for ( i in 1:10000){
  v = X %*% b
  # fit spline model
  spfit <- smooth.spline(x=v, y=y)
  g <- predict(spfit, v, deriv=0)$y
  gprime <- predict(spfit, v, deriv=1)$y
  
  # Update weight
  w <- diag(as.vector(gprime)/sum(gprime^2))
  bhat <- v + (y - g)/gprime
  
  # Update beta with weighted least squares
  b_new <- solve(t(X) %*% w %*%X) %*% t(X) %*% w %*% as.vector(bhat)
  
  if (sqrt(sum((b-b_new)^2))<1e-6){
    print(b_new)
    break
  } else {
    b <- b_new
    print(i)
    next
  }
  
}

## Compare with ppr function results
b_new
pprfit$alpha


