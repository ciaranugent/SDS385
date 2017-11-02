wi <- function(z){1/(1+exp(-z))}

# generate data

p <- 6
n <- 1000

beta <- as.matrix(sample(seq(-7,7,.2),p),p,1)
W <- rep(1,n)
X <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1,byrow=TRUE))
y <- rbinom(n,W,wi(crossprod(t(X),beta)))

# set step size

alpha <- c(.1,0.01,0.001,0.0001)

# Basic SGD

all.betas <- NULL
ls        <- NULL
for (h in 1:length(alpha)){
	betas <- NULL
	l <- NULL
	beta0 <- as.matrix(rep(0,p),p,1)
	for (i in 1:10000){
		for (k in 1:10){
			j <- sample(1:n,1)
			beta0 <- beta0 - alpha[h]*W[j]*(wi(sum(X[j,]*beta0))-y[j])*X[j,]
		}
	betas <- rbind(betas,t(beta0))
	l <- c(l)
	}
	all.betas[[h]] <- betas
}

# for my run I had a true beta of 
# beta
# -0.8
#  6.6
# -0.2
#  3.4
#  6.8
# -6.4

b.num <- 1

par(mfrow=c(2,2))
plot(seq(1,100000,10),all.betas[[1]][,b.num],main='alpha = 0.1',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[1]][,b.num]),max(all.betas[[1]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,100000,10),all.betas[[2]][,b.num],main='alpha = 0.01',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[1]][,b.num]),max(all.betas[[1]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,100000,10),all.betas[[3]][,b.num],main='alpha = 0.001',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[1]][,b.num]),max(all.betas[[1]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,100000,10),all.betas[[4]][,b.num],main='alpha = 0.0001',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[1]][,b.num]),max(all.betas[[1]][,b.num])))
abline(beta[b.num],0,col=2)

# Write a function for the exponentially weighted moving average

mov.avg <- function(x,y,lambda){
	lambda*x + (1-lambda)*y
}

# Write a function for decaying step size

st.size <- function(C,t,t0,a){C*(t+t0)^(-a)}
alpha <- c(0.6,0.6,0.85,0.85)
Cs <- c(1,10,1,10)

all.betas <- NULL
ls        <- NULL
for (h in 1:length(alpha)){
	betas <- NULL
	l <- NULL
	beta0 <- as.matrix(rep(0,p),p,1)
	for (i in 1:20000){
		for (k in 1:10){
			j <- sample(1:n,1)
			beta0 <- beta0 - st.size(Cs[h],(i-1)*10 + k,1,alpha[h])*W[j]*(wi(sum(X[j,]*beta0))-y[j])*X[j,]
		}
	betas <- rbind(betas,t(beta0))
	l <- c(l)
	}
	all.betas[[h]] <- betas
}

b.num <- 1

par(mfrow=c(2,2))
plot(seq(1,200000,10),all.betas[[1]][,b.num],main='alpha = 0.6, C = 1',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[2]][,b.num]),max(all.betas[[2]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,200000,10),all.betas[[2]][,b.num],main='alpha = 0.6, C = 10',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[2]][,b.num]),max(all.betas[[2]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,200000,10),all.betas[[3]][,b.num],main='alpha = 0.85, C = 1',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[2]][,b.num]),max(all.betas[[2]][,b.num])))
abline(beta[b.num],0,col=2)

plot(seq(1,200000,10),all.betas[[4]][,b.num],main='alpha = 0.85, C = 10',xlab='Iteration',ylab='Value',ylim=c(min(all.betas[[2]][,b.num]),max(all.betas[[2]][,b.num])))
abline(beta[b.num],0,col=2)


