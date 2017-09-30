library(microbenchmark)

# Methods

inverse.method  <- function(a,b){solve(a)%*%b}

default.method <- function(a,b){solve(a,b)}

cholesky.method <- function(a,b){temp <- chol(a)
								    Z <- backsolve(temp,b,transpose=TRUE)
								    backsolve(temp,Z)}

# LU.method       <- function(a,b){ex <- lu.decomposition(a)
								 # Z  <- forwardsolve(ex$L,b)
								 # backsolve(ex$U,Z)}

n <- 5000

p.seq <- seq(100,1400,100)
results <- matrix(0,length(p.seq),9)

for(i in 1:length(p.seq)){
p <- p.seq[i]

W <- rep(1,n)
X <- matrix(rnorm(n*p),n,p,byrow=TRUE)
y <- matrix(rnorm(n),n,1)

A <- t(X)%*%(W*X)
B <- t(X)%*%(W*y)

temp <- microbenchmark(inverse.method(A,B),cholesky.method(A,B),default.method(A,B),times=1000)

results[i,1] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,2] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,3] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
results[i,4] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,5] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,6] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
results[i,7] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,8] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,9] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
}

# Results
             # [,1]        [,2]        [,3]       [,4]        [,5]        [,6]        [,7]       [,8]        [,9]
 # [1,]    21839.53    22927.01    31932.65    1907.51    2085.855    3723.448    4795.618    5060.63    8120.759
 # [2,]   165835.96   177967.96   209017.39   12614.88   14915.455   23108.213   28016.838   31078.40   45345.108
 # [3,]   553638.66   590374.05   648068.58   41783.80   46237.540   64488.395   82185.391   89354.42  115522.362
 # [4,]  1325783.71  1384546.81  1462069.67   99249.99  108410.885  138230.378  178955.770  196592.70  230566.798
 # [5,]  2615162.53  2690627.95  2790617.87  197151.81  217704.590  254854.806  333649.173  366367.84  416999.629
 # [6,]  4330101.09  4622264.01  4820267.34  343707.22  366656.940  418873.628  569034.237  613013.40  680169.400
 # [7,]  6966484.40  7257330.57  7686232.63  545405.04  578434.660  646068.377  885818.876  948818.60 1021479.847
 # [8,] 10818630.44 11037032.14 11367222.05  822436.40  857049.210  927941.556 1292864.145 1380397.69 1491179.251
 # [9,] 15517891.67 15815401.93 16215738.65 1174276.62 1216939.525 1306508.492 1815652.826 1931693.08 2072103.061
# [10,] 21351285.30 21735333.99 22207324.03 1617636.30 1681361.505 1785838.446 2468993.664 2605513.96 2765471.715
# [11,] 28447447.22 28904397.33 29458970.03 2156167.05 2231135.900 2351091.428 3266444.886 3423474.92 3627253.897
# [12,] 36920104.09 37462812.58 38044814.17 2824675.47 2908935.950 3018366.681 4209451.067 4387706.99 4619942.334
# [13,] 46934690.47 47531433.48 48160762.71 3593203.86 3693188.065 3836789.305 5334147.692 5543569.88 5783250.116
# [14,] 58599303.95 59390374.02 62404234.11 4517585.76 4648198.395 4864046.863 6648286.035 6906759.67 7256140.249

plot(p.seq[1:14],results[1:14,2],type='l',main='N = 5000 Observations',xlab='Number of P Variables',ylab='Time (microseconds)')
lines(p.seq[1:14],results[1:14,1],type='l',lty=2)
lines(p.seq[1:14],results[1:14,3],type='l',lty=2)
lines(p.seq[1:14],results[1:14,5],type='l',lty=1,col=2)
lines(p.seq[1:14],results[1:14,4],type='l',lty=2,col=2)
lines(p.seq[1:14],results[1:14,6],type='l',lty=2,col=2)
lines(p.seq[1:14],results[1:14,8],type='l',lty=1,col=3)
lines(p.seq[1:14],results[1:14,7],type='l',lty=2,col=3)
lines(p.seq[1:14],results[1:14,9],type='l',lty=2,col=3)

n <- 1000

p.seq <- seq(100,500,50)
results <- matrix(0,length(p.seq),9)

for(i in 1:length(p.seq)){
p <- p.seq[i]

W <- rep(1,n)
X <- matrix(rnorm(n*p),n,p,byrow=TRUE)
y <- matrix(rnorm(n),n,1)

A <- crossprod(X,W*X)
B <- crossprod(X,W*y)

temp <- microbenchmark(inverse.method(A,B),cholesky.method(A,B),default.method(A,B),times=1000)

results[i,1] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,2] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,3] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
results[i,4] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,5] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,6] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
results[i,7] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[1]/1000
results[i,8] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[2]/1000
results[i,9] <- quantile(temp[[2]][temp[[1]]=='default.method(A, B)'],probs=c(.05,.5,.95))[3]/1000
}

plot(p.seq,results[,2],type='l',main='N = 1000 Observations',xlab='Number of P Variables',ylab='Time (microseconds)')
lines(p.seq,results[,1],type='l',lty=2)
lines(p.seq,results[,3],type='l',lty=2)
lines(p.seq,results[,5],type='l',lty=1,col=2)
lines(p.seq,results[,4],type='l',lty=2,col=2)
lines(p.seq,results[,6],type='l',lty=2,col=2)
lines(p.seq,results[,8],type='l',lty=1,col=3)
lines(p.seq,results[,7],type='l',lty=2,col=3)
lines(p.seq,results[,9],type='l',lty=2,col=3)
legend(125,250000,c('Inversion Method','Cholesky Decomposition','Default R'),lty=1,col=c(1:3))

#############################################
############## Sparse Matrices ##############
#############################################

library(Matrix)

# Make a Sparse Matrix

N <- 2000

y <- matrix(rnorm(N),N,1)

p.seq <- seq(100,1000,100)

results <- matrix(0,length(p.seq),18)

for (i in 1:length(p.seq)){
dat      <- rnorm(N*p.seq[i])
mask     <- rbinom(N*p.seq[i],1,0.05)
X.reg    <- matrix(dat*mask, nrow=N)
X.sparse <- Matrix(dat*mask, nrow=N)

A.reg    <- crossprod(X.reg,X.reg) 
B.reg    <- crossprod(X.reg,y)

A.sparse <- crossprod(X.sparse,X.sparse) 
B.sparse <- crossprod(X.sparse,y)
	
temp <- microbenchmark(inverse.method(A.reg,B.reg),cholesky.method(A.reg,B.reg),default.method(A.reg,B.reg),inverse.method(A.sparse,B.sparse),cholesky.method(A.sparse,B.sparse),default.method(A.sparse,B.sparse),times=1000)

results[i,1] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[1]/1000
results[i,2] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[2]/1000
results[i,3] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[3]/1000
results[i,4] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[1]/1000
results[i,5] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[2]/1000
results[i,6] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[3]/1000
results[i,7] <- quantile(temp[[2]][temp[[1]]=='default.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[1]/1000
results[i,8] <- quantile(temp[[2]][temp[[1]]=='default.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[2]/1000
results[i,9] <- quantile(temp[[2]][temp[[1]]=='default.method(A.reg, B.reg)'],probs=c(.05,.5,.95))[3]/1000

results[i,10] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[1]/1000
results[i,11] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[2]/1000
results[i,12] <- quantile(temp[[2]][temp[[1]]=='inverse.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[3]/1000
results[i,13] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[1]/1000
results[i,14] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[2]/1000
results[i,15] <- quantile(temp[[2]][temp[[1]]=='cholesky.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[3]/1000
results[i,16] <- quantile(temp[[2]][temp[[1]]=='default.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[1]/1000
results[i,17] <- quantile(temp[[2]][temp[[1]]=='default.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[2]/1000
results[i,18] <- quantile(temp[[2]][temp[[1]]=='default.method(A.sparse, B.sparse)'],probs=c(.05,.5,.95))[3]/1000
}

plot(p.seq,results[,2],type='l',main='N = 2000 Observations and Sparsity = 0.05',xlab='Number of P Variables',ylab='Time (microseconds)')
lines(p.seq,results[,1],type='l',lty=2)
lines(p.seq,results[,3],type='l',lty=2)
lines(p.seq,results[,5],type='l',lty=1,col=2)
lines(p.seq,results[,4],type='l',lty=2,col=2)
lines(p.seq,results[,6],type='l',lty=2,col=2)
lines(p.seq,results[,8],type='l',lty=1,col=3)
lines(p.seq,results[,7],type='l',lty=2,col=3)
lines(p.seq,results[,9],type='l',lty=2,col=3)
lines(p.seq,results[,11],type='l',col=4)
lines(p.seq,results[,10],type='l',lty=2,col=4)
lines(p.seq,results[,12],type='l',lty=2,col=4)
lines(p.seq,results[,14],type='l',lty=1,col=5)
lines(p.seq,results[,13],type='l',lty=2,col=5)
lines(p.seq,results[,15],type='l',lty=2,col=5)
lines(p.seq,results[,17],type='l',lty=1,col=6)
lines(p.seq,results[,16],type='l',lty=2,col=6)
lines(p.seq,results[,18],type='l',lty=2,col=6)
legend(125,2500000,c('Inversion Method','Cholesky Decomposition','Default R','Inversion Method  - Sparse','Cholesky Decomposition - Sparse','Default R - Sparse'),lty=1,col=c(1:6))


            # [,1]        [,2]        [,3]        [,4]        [,5]        [,6]        [,7]       [,8]        [,9]       [,10]       [,11]      [,12]
 # [1,]    2206.774    2348.119    3407.088    200.5694    234.5615    465.2911    492.3897    535.303    943.9588    5822.209    7117.582   11659.51
 # [2,]   16752.454   18759.775   24114.837   1278.1715   1534.7835   3053.3930   2830.6324   3235.744   5421.2931   20242.844   24562.355   42360.52
 # [3,]   55154.220   60339.024   90329.198   4171.1936   4702.2100   7702.5932   8215.1421   9118.687  15290.0273   48633.698   54007.957  112514.97
 # [4,]  148759.082  161104.903  246331.589  10944.5916  12827.6640  18954.8687  20051.4051  22901.587  31148.8596  119554.355  132333.193  214285.38
 # [5,]  270009.917  300012.370  385986.860  20100.3356  22998.3895  29422.1671  35138.7919  40155.516  51447.2349  189460.013  222188.160  381495.53
 # [6,]  471223.032  542366.490  661833.831  35604.2455  40615.7025  49232.2959  61335.3596  69125.503  80717.0220  313526.957  384613.628  532682.53
 # [7,]  740649.438  855357.816 1248341.649  56857.2733  62851.6525  83068.6097  92980.3924 103097.145 145404.9545  472651.721  583874.919  962244.10
 # [8,] 1110944.773 1237826.046 1744901.627  86078.1756  95243.5850 124374.0806 135250.1553 151566.223 213003.2560  696287.038  807815.525 1350860.54
 # [9,] 1630823.702 1864247.722 2440534.759 123598.7984 134329.9440 159086.1372 191686.5862 213039.993 268860.0081  989164.514 1250947.390 1709369.05
# [10,] 2198006.974 2514805.979 3081209.370 168287.8706 183609.9510 216176.3373 259414.3986 285231.070 356146.7791 1355666.421 1647220.569 2248802.89
           # [,13]      [,14]     [,15]      [,16]     [,17]     [,18]
 # [1,]   269.8023   374.4225  1157.418   37.11630   58.7260  148.1164
 # [2,]   619.6907  1230.0645  3764.521   69.03455  136.3105  423.9732
 # [3,]  1111.7552  2130.3790  6095.793  127.62600  211.9045  664.6326
 # [4,]  2357.9482  4795.3015 11488.578  266.32710  385.6535 1138.0525
 # [5,]  3383.2840  6182.0315 13042.042  372.49285  468.0460 1293.3976
 # [6,]  5370.5189  9242.8215 17674.086  521.02735  616.2310 1767.0629
 # [7,]  7239.1132 12047.4585 24846.561  679.94410  771.3720 2222.6820
 # [8,]  8988.0461 15416.0625 32493.512  855.46415  948.5350 3205.4267
 # [9,] 12793.6726 19575.9695 37346.474 1044.78170 1212.7105 3766.2317
# [10,] 16478.2961 24967.1440 42896.800 1250.46200 1429.4510 3809.7956