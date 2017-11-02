library(glmnet)

diabY <- as.matrix(read.csv('https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesY.csv',header=FALSE))
diabX <- as.matrix(read.csv('https://raw.githubusercontent.com/jgscott/SDS385/master/data/diabetesX.csv'))

temp <- glmnet(diabX,diabY,"gaussian")
plot(temp,xvar='lambda')

path <- temp$lambda
vals <- colSums((matrix(rep(diabY,100),byrow=FALSE,dim(diabY)[1],100) - predict(temp,diabX,s=c(path)))^2)/dim(diabY)[1]

plot(log(path),vals,type='l',xlab='log lambda',ylab='MSE',main='LASSO Fit')

####################################
######### Cross Validation #########
####################################

# N-fold Cross Validation
# Inputs: data, sample size, number of folds, bandwith, Kernal

cv.N <- function(n,N,X,Y,lams){	
	pred.err <- NULL
	
	draws <- NULL
	for (i in 1:N){
		draws <- c(draws,rep(i,floor(n/N)))
	}
	draws <- c(draws,sample(1:N,n-N*floor(n/N)))
	samp  <- sample(draws,n)
	
	for (s in 1:N){
		ytrain <- as.matrix(Y[-which(samp==s),])
		ytest  <- as.matrix(Y[which(samp==s),])
		xtrain <- X[-which(samp==s),]
		xtest  <- X[which(samp==s),]

		temp2 <- glmnet(xtrain,ytrain,"gaussian")		
		pred.err <- cbind(pred.err,colSums((matrix(rep(ytest,length(lams)),byrow=FALSE,dim(ytest)[1],length(lams)) - predict(temp2,xtest,s=c(lams)))^2)/dim(ytest)[1])
	}
	pred.err
}


alphas <- seq(0,50,.1)

MSES <- cv.N(dim(diabY)[1],5,diabX,diabY,alphas)

MSES <- cv.N(dim(diabY)[1],5,diabX,diabY,alphas)
AvgMSE <- rowMeans(MSES)
plot(log(alphas),AvgMSE,type='l',main="Average Test Error - 5-fold Cross-Validation",xlab="Log Lambda",ylab="Test Error")


####################################
######### Mallow's Cp #########
####################################

lmmod <- lm(diabY ~ diabX)
sig   <- sum((diabY - predict(lmmod,as.data.frame(diabX)))^2)/(dim(diabX)[1]-dim(diabX)[2]-1)

Cp <- colSums((matrix(rep(diabY,100),byrow=FALSE,dim(diabY)[1],100) - predict(temp,diabX,s=c(path)))^2)/dim(diabY)[1] + 2*temp$df*sig/dim(diabX)[1]
plot(log(path),Cp,type='l',main="Mallow's Cp",xlab="Log Lambda",ylab="Cp Statistic",col=3)


# Plot All

plot(log(path),vals,type='l',xlab='log lambda',ylab='Error',main='Comparison')
points(log(path),Cp,type='l',col=3)
points(log(alphas),AvgMSE,type='l',col=2)

legend(-4.5,5500,c("LASSO Fit","5-Fold Cross-Validation","Mallow's Cp"),lty=1,col=1:3)