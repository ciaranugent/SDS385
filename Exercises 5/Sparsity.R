########################
######## Part A ########
########################

Sy <- function(lambda,y){ifelse((abs(y)-lambda) < 0,0,sign(y)*(abs(y)-lambda))}
Hy <- function(lambda,y){ifelse(abs(y) > lambda,y,0)}

vals <- seq(-5,5,.01)

par(mfrow=c(2,2))
plot(vals,Sy(2,vals),main='Soft Thresholding: lambda = 2',xlab='y',ylab='S(y)')
plot(vals,Hy(2,vals),main='Hard Thresholding: lambda = 2',xlab='y',ylab='H(y)')
plot(vals,Sy(1,vals),main='Soft Thresholding: lambda = 1',xlab='y',ylab='S(y)')
plot(vals,Hy(1,vals),main='Hard Thresholding: lambda = 1',xlab='y',ylab='H(y)')

########################
######## Part B ########
########################

N <- 100
j <- 0.25

thetas <- rnorm(N,0,5)*rbinom(N,1,j)
sigmas <- rgamma(N,1)

zs <- rnorm(N,thetas,sqrt(sigmas))

plot(thetas,Sy(2,zs),ylab='Soft Thresold Estimate of Theta',xlab='True Theta',main='Soft Thresholding: 75% Sparse',ylim=c(-10,10))
points(thetas,Sy(1,zs),col=2)
points(thetas,Sy(.5,zs),col=3)
points(thetas,Sy(.1,zs),col=4)
abline(0,1)
legend(5,-4,c('lambda = 2','lambda = 1','lambda = 0.5','lambda = 0.1'),pch=1,col=1:4)

lambdas <- seq(0,5,.1)
MSES     <- NULL

sparsity <- c(0.5,0.4,0.3,0.2)

for (j in 1:length(sparsity)){
	thetas <- rnorm(N,0,5)*rbinom(N,1,sparsity[j])
	zs <- rnorm(N,thetas,sqrt(sigmas))<-
	MSE     <- NULL
	for (i in 1:length(lambdas)){
		MSE <- c(MSE,mean((thetas-Sy(lambdas[i],zs))^2)) 
	}
	MSES <- rbind(MSES,MSE)
}

plot(lambdas,MSES[1,],type='l',ylim=c(min(MSES),max(MSES)),xlab='Lambda',ylab='MSE',main='Plots of MSE for Varying Levels of Sparsity')

for(i in 2:4){
	points(lambdas,MSES[i,],type='l',col=i)
}

legend(0,6,c('0.5','0.6','0.7','0.8'),lty=1,col=1:4)

