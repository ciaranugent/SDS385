my.data<-read.csv('https://raw.githubusercontent.com/jgscott/SDS385/master/data/wdbc.csv')

X <- as.matrix(cbind(1,my.data[,3:12]))
y <- ifelse(my.data[,2]=='M',1,0)

wi <- function(z){1/(1+exp(-z))}

# Initialize beta
beta <- matrix(rep(0,11),11,1)
betas <- NULL

for (i in 1:10){
	z <- X%*%beta
	beta <- beta - solve(crossprod(X,diag(c(wi(z)*(1-wi(z))))%*%X),crossprod(X,wi(z)-y))
	betas <- cbind(betas,beta)
}


# Solve Weighted Least Squares form from part C, I have added in a jitter of 0.0000001 to prevent dividing by zero, so the answers are the same as above upto a certain number of decimal places
for (i in 1:10){
	z <- X%*%beta
	beta <- solve(crossprod(X,crossprod(-diag(c(wi(z)*(1-wi(z)))),X)),crossprod(X,-diag(c(wi(z)*(1-wi(z)))))%*%(-1/(1-wi(z)+0.0000001)+(1/(wi(z)*(1-wi(z))+0.0000001))*y + z))
	betas <- cbind(betas,beta)
}

