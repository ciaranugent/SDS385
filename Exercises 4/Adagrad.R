library(Matrix)
library(Rcpp)
library(RcppEigen)

wi <- function(z){1/(1+exp(-z))}

# generate data

p <- 6
n <- 1000
jitter <- 0.000001

beta <- as.matrix(sample(seq(-7,7,.2),p),p,1)
W <- rep(1,n)
temp <- rnorm(n*p)*rbinom(n*p,1,0.05)
X <- matrix(temp,p,n,byrow=FALSE)
y <- rbinom(n,1,as.numeric(wi(crossprod(X,beta))))

X <- Matrix(temp,n,p,byrow=TRUE)

# Initialize Values

n <- dim(X)[1]
alpha <- 0.01
sm.factor <- 0.00000001
h_grad <- 0
beta0 <- as.matrix(rep(0,p),p,1)
lambda <- 0.0001

# choose only rows that have non-zero values

vals <- unique(summary(X)[,1])

# Adagrad

betas <- NULL
l <- NULL
h_grad <- 0
beta0 <- as.matrix(rep(0,p),p,1)
for (i in 1:10000){
	for (k in 1:10){
		j <- sample(vals,1)
		grad <- W[j]*(wi(sum(X[j,]*beta0))-y[j])*X[j,] + lambda*sign(beta0)
		h_grad <- h_grad + grad^2
		a_grad <- grad/(jitter+sqrt(h_grad))
		beta0 <- beta0 - alpha*a_grad
	}
betas <- rbind(betas,t(beta0))
}






