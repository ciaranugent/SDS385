my.data<-read.csv('https://raw.githubusercontent.com/jgscott/SDS385/master/data/wdbc.csv')

X <- scale(as.matrix(my.data[,3:12]))
X <- cbind(1,X)
y <- ifelse(my.data[,2]=='M',1,0)

wi <- function(z){1/(1+exp(-z))}

betas <- NULL
beta0 <- as.matrix(rnorm(dim(X)[2]),dim(X)[2],1)

jit <- 0.000000001

par(mfrow=c(2,2))


ls <- NULL
alpha <- .05
for (i in 1:20000){
	beta0 <- beta0 - alpha*crossprod(X,wi(crossprod(t(X),beta0))-y)
	betas <- rbind(betas,t(beta0))
	l <- -sum(y*log(jit+wi(crossprod(t(X),beta0))))-sum((1-y)*log(jit + 1-wi(crossprod(t(X),beta0))))
	ls <- c(ls,as.numeric(l))
}

plot(ls,ylab='Negative Log-Likelihood',xlab='Iteration',main='alpha=0.05',ylim=c(0,350))


beta0 <- as.matrix(rnorm(dim(X)[2]),dim(X)[2],1)
ls <- NULL
alpha <- .01
for (i in 1:20000){
	beta0 <- beta0 - alpha*crossprod(X,wi(crossprod(t(X),beta0))-y)
	betas <- rbind(betas,t(beta0))
	l <- -sum(y*log(jit+wi(crossprod(t(X),beta0))))-sum((1-y)*log(jit + 1-wi(crossprod(t(X),beta0))))
	ls <- c(ls,as.numeric(l))
}

plot(ls,ylab='Negative Log-Likelihood',xlab='Iteration',main='alpha=0.01',ylim=c(0,350))

beta0 <- as.matrix(rnorm(dim(X)[2]),dim(X)[2],1)
ls <- NULL
alpha <- .005
for (i in 1:20000){
	beta0 <- beta0 - alpha*crossprod(X,wi(crossprod(t(X),beta0))-y)
	betas <- rbind(betas,t(beta0))
	l <- -sum(y*log(jit+wi(crossprod(t(X),beta0))))-sum((1-y)*log(jit + 1-wi(crossprod(t(X),beta0))))
	ls <- c(ls,as.numeric(l))
}

plot(ls,ylab='Negative Log-Likelihood',xlab='Iteration',main='alpha=0.005',ylim=c(0,350))

beta0 <- as.matrix(rnorm(dim(X)[2]),dim(X)[2],1)
ls <- NULL
alpha <- .001
for (i in 1:20000){
	beta0 <- beta0 - alpha*crossprod(X,wi(crossprod(t(X),beta0))-y)
	betas <- rbind(betas,t(beta0))
	l <- -sum(y*log(jit+wi(crossprod(t(X),beta0))))-sum((1-y)*log(jit + 1-wi(crossprod(t(X),beta0))))
	ls <- c(ls,as.numeric(l))
}

plot(ls,ylab='Negative Log-Likelihood',xlab='Iteration',main='alpha=0.005',ylim=c(0,350))

