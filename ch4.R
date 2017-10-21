###  R code for Chapter 4  
require(dlm)


###
### Section 1: MLE
###


### MLE -- Lake Superior data

y <- ts(read.table("Datasets/lakeSuperior.dat",
                   skip=3)[, 2], start = c(1900, 1))
build <- function(parm) {
  dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
}
fit <- dlmMLE(y, rep(0,2), build)
fit$convergence
unlist(build(fit$par)[c("V","W")])
## ask for Hessian
fit <- dlmMLE(y, rep(0,2), build, hessian=T)
avarLog <- solve(fit$hessian)
avar <- diag(exp(fit$par)) %*% avarLog %*%
  diag(exp(fit$par)) # Delta method
sqrt(diag(avar)) # estimated standard errors
require(nlme)
avar1 <- solve(fdHess(exp(fit$par), function(x)
                      dlmLL(y,build(log(x))))$Hessian)
sqrt(diag(avar1))
## direct parametrization
build <- function(parm) {
  dlmModPoly(order=1, dV=parm[1], dW=parm[2])
}
fit <- dlmMLE(y, rep(0.23,2), build, lower=c(1e-6,0), hessian=T)
fit$convergence
unlist(build(fit$par)[c("V","W")])
avar <- solve(fit$hessian)
sqrt(diag(avar))

StructTS(y,"level")


### The following example refers to Chapter 3, page 130:
### Denmark and Spain investments -- integrated random walk SUTSE

invest <- ts(matrix(scan("Datasets/invest2.dat"), nc = 2, byrow = TRUE),
             start = 1960, names = c("Denmark", "Spain"))
invest[, 2] <- invest[, 2] / 100
build <- function(par) {
    m <- dlmModPoly(2)
    m$FF <- m$FF %x% diag(2)
    m$GG <- m$GG %x% diag(2)
    W1 <- diag(1e-7, nr = 2)
    W2 <- diag(exp(par[1 : 2]))
    W2[1,2] <- W2[2,1] <- tanh(par[3]) *
        prod(exp(0.5 * par[1 : 2]))
    m$W <- bdiag(W1, W2)
    V <- diag(exp(par[4 : 5]))
    V[1,2] <- V[2,1] <- tanh(par[6]) *
        prod(exp(0.5 * par[4 : 5]))
    m$V <- V
    m$m0 <- rep(0,4)
    m$C0 <- diag(4) * 1e7
    return(m)
}

mle <- vector("list", 16)
init <- cbind(as.matrix(expand.grid(rep(list(c(-3, 3)), 4))),
              matrix(0, nc = 2, nr = 16))
init <- init[,c(1, 2, 5, 3, 4, 6)]; dimnames(init) <- NULL
for (i in 1:length(mle)) {
    mle[[i]] <- try(dlmMLE(invest, parm = init[i,], build = build,
                           control = list(trace = 3, REPORT = 2)))
}
good <- seq(along=mle)[sapply(mle, function(x) length(x) > 1 && x$conv == 0)]
mleX <- sapply(mle[good], function(x) unlist(x[1 : 2]))
mleX <- mleX[, order(mleX["value", ])]

fit <- build(mleX[-NROW(mleX), 1])
### etc...



###
###  Section 2: Bayesian inference
###

### Section 4.3.2: Conjugate inference with discount factors

# require function "dlmFilterDF"

y <- ts(read.table("Datasets/lakeSuperior.dat",
                   skip = 3)[, 2], start = c(1900, 1))
mod <- dlmModPoly(1, dV = 1)
modFilt <- dlmFilterDF(y, mod, DF = 0.9)
beta0 <- 20; alpha0 <- 2

par(mfrow = c(3, 1), mar = c(1.5, 4, 0.5, 0.5) + 0.1, cex = 0.6)
## Filtering estimates
out <- residuals(modFilt)
beta <- beta0 + cumsum(out$res^2) / 2
alpha <- alpha0 + (1 : length(y)) / 2
Ctilde <- unlist(dlmSvd2var(modFilt$U.C, modFilt$D.C))[-1]
prob <- 0.95
tt <- qt(prob, df = 2 * alpha)
lower <-  dropFirst(modFilt$m) - tt * sqrt(Ctilde * beta / alpha)
upper <-  dropFirst(modFilt$m) + tt * sqrt(Ctilde * beta / alpha)
plot(y, ylab = "Filtering level estimates", type = "o", 
		ylim = c(18, 40), col = "darkgray")
lines(dropFirst(modFilt$m), type = "o")
lines(lower, lty = 2, lwd = 2)
lines(upper, lty = 2, lwd = 2)
## One-step-ahead forecasts
sigma2 <- c(beta0 / (alpha0-1), beta / (alpha-1))
Qt <- out$sd^2 * sigma2[-length(sigma2)]
alpha0T = c(alpha0,alpha)
tt <- qt(prob, df = 2 * alpha0T[-length(alpha0T)])
parf <- c(beta0 / alpha0, beta / alpha)
parf <- parf[-length(parf)] * out$sd^2
lower <- dropFirst(modFilt$f) - tt * sqrt(parf)
upper <- dropFirst(modFilt$f) + tt * sqrt(parf)
plot(y, ylab = "One-step-ahead forecasts", type = "o",
     ylim = c(20, 40), col = "darkgray")
lines(window(modFilt$f, start=1902),  type="o") 
lines(lower, lty = 2, lwd = 2)
lines(upper, lty = 2, lwd = 2)
## Smoothing estimates
modFilt$mod$JW <- matrix(1)
X <- unlist(dlmSvd2var(modFilt$U.W, modFilt$D.W))[-1]
modFilt$mod$X <- matrix(X)
modSmooth <- dlmSmooth(modFilt)
Stildelist <- dlmSvd2var(modSmooth$U.S, modSmooth$D.S)
TT <- length(y)
pars <- unlist(Stildelist) * (beta[TT] / alpha[TT])
tt <- qt(prob, df = 2 * alpha[TT])
plot(y, ylab = "Smoothing level estimates", 
     type = "o", ylim = c(20, 40), col = "darkgray")
lines(dropFirst(modSmooth$s),  type = "o")
lines(dropFirst(modSmooth$s - tt * sqrt(pars)), lty = 3, lwd = 2)
lines(dropFirst(modSmooth$s + tt * sqrt(pars)), lty = 3, lwd = 2)

## Choosing the discount factor
vecDF <- c(1,0.9,0.8,0.3)
mat <- matrix(0, nrow=length(vecDF), ncol=5,dimnames = 
	list(c("", "","",""), c("DF", "MAPE", "MAD","MSE","SD")))
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(2,2)) 
for (i in 1:length(vecDF)){
   mod <- dlmModPoly(1,dV=1)
   modFilt <- dlmFilterDF(y, mod, DF=vecDF[i])
   out <- residuals(modFilt)
   beta0 <- 20
   alpha0 <-2
   beta <- beta0 + 1/2 * cumsum(out$res^2)
   alpha <- alpha0 + (1: length(y))*(1/2)
   sigmaqhat <- beta/(alpha-1) 
   sigmaqhat <- beta/(alpha-1) 
   # Plot the estimates of sigma^2
   # ts.plot(sigmaqhat, main=paste("DF=", as.character(vecDF[i])), ylab="")  
   # Plot the one-step-ahead forecasts   
   plot(window(y, start=1960), ylab="", ylim=c(20,45), lty=1, 
		main=paste("DF=", as.character(vecDF[i])), col="darkgray") 
   lines(window(modFilt$f, start=1960), lty=1)  
   mat[i,1] <- vecDF[i]
   mat[i,2] <- mean(abs(modFilt$f - y)/y)
   mat[i,3] <- mean(abs(modFilt$f - y))
   mat[i,4] <- mean((modFilt$f - y)^2)
   mat[i,5] <- sigmaqhat[length(y)]
}

round(mat,4)


###  Section 4.3.3. 

y <- ts(read.table("Datasets/lakeSuperior.dat",
                   skip = 3)[, 2], start = c(1900, 1))
beta0 <- 20; alpha0 <- 2 ; TT <- length(y)
mod <- dlmModPoly(1, dV = 1)
modFilt <- dlmFilterDF(y, mod, DF = 0.9)
out <- residuals(modFilt)   
DFstar <- 0.9
delta <- DFstar^(0 : TT)
alpha <- delta[-1] * alpha0 + cumsum(delta[-TT]) / 2
res <- as.vector(out$res) 
beta <- delta[-1] * beta0
for (i in 1 : TT)
    beta[i] <- beta[i] + 0.5 * sum(delta[i:1] * res[1:i]^2)
alphaStar <- DFstar * c(alpha0, alpha)
betaStar <- DFstar * c(beta0, beta)
tt <- qt(0.95, df = 2 * alphaStar[1 : TT])
param <- sqrt(out$sd) * (betaStar / alphaStar)[1 : TT]

par(mar = c(1.5, 4, 0.5, 0.5) + 0.1, cex = 0.6)
plot(y, ylab = "Observed/One step ahead forecasts",
     type = "o", col = "darkgray", ylim =c(20, 45)) 
lines(window(modFilt$f, start = 1902), type = "o")
lines(window(modFilt$f, start = 1902) - tt * sqrt(param),
      lty = 2, lwd = 2)
lines(window(modFilt$f, start = 1902) + tt * sqrt(param),
      lty = 2, lwd = 2)



###
### Section 3: Simulation-based Bayesian inference
###


### Illustration: Gibbs for local level + Nile

a1 <- 2 
b1 <- 0.0001
a2 <- 2
b2 <- 0.0001
## starting values
psi1 <- 1
psi2 <- 1
mod_level <- dlmModPoly(1, dV = 1 / psi1, dW = 1 / psi2)

mc <- 1500
psi1_save <- numeric(mc)
psi2_save <- numeric(mc)
n <- length(Nile)
sh1 <- a1 + n / 2
sh2 <- a2 + n / 2
set.seed(10)

for (it in 1 : mc)
{
    ## draw the states: FFBS
    filt <- dlmFilter(Nile, mod_level)
    level <- dlmBSample(filt)
    ## draw observation precision psi1
    rate <- b1 + crossprod(Nile - level[-1]) / 2
    psi1 <- rgamma(1, shape = sh1, rate = rate)
    ## draw system precision psi2
    rate <- b2 + crossprod(level[-1] - level[-n]) / 2
    psi2 <- rgamma(1, shape = sh2, rate = rate)
    ## update and save
    V(mod_level) <- 1 / psi1
    W(mod_level) <- 1 / psi2
    psi1_save[it] <- psi1
    psi2_save[it] <- psi2
}


### Section 4.5.1: Gibbs sampling for "d-inverse-gamma" model

invSpain <- ts(read.table("Datasets/invest2.dat",
                          colClasses = "numeric")[,2]/1000, start = 1960)
set.seed(5672)
MCMC <- 12000
gibbsOut <- dlmGibbsDIG(invSpain, mod = dlmModPoly(2), a.y = 1, b.y = 1000,
                        a.theta = 10, b.theta = 1000, n.sample = MCMC,
                        thin = 1, save.states = FALSE)

## output analysis
burn <- 2000

## Convergence assessment and autocorrelations
par(mfrow=c(3,3), mar=c(3.1,2.1,2.1,1.1))
plot(gibbsOut$dV[-(1:burn)], cex=0.5, xlab="", ylab="", main=expression(V))
plot(gibbsOut$dW[-(1:burn),1], cex=0.5, xlab="", ylab="", main=expression(W[11]))
plot(gibbsOut$dW[-(1:burn),2], cex=0.5, xlab="", ylab="", main=expression(W[22]))
use <- MCMC - burn
from <- 0.05 * use
plot(ergMean(gibbsOut$dV[-(1:burn)], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$dW[-(1:burn),1], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$dW[-(1:burn),2], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
acf(gibbsOut$dV[-(1:burn)], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$dW[-(1:burn),1], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$dW[-(1:burn),2], ci=0, ylim=c(0,1), ylab="", main="")

## posterior estimates of variances, with estimated MC standard error
mcmcMeans(cbind(gibbsOut$dV[-(1:burn),], gibbsOut$dW[-(1:burn),]))

par(mfrow=c(1,3), mar=c(5.1,4.1,1,1))
plot(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),1],
     xlab=expression(V), ylab=expression(W[11]), pch='.')
plot(gibbsOut$dV[-(1:burn)], gibbsOut$dW[-(1:burn),2],
     xlab=expression(V), ylab=expression(W[22]), pch='.')
plot(gibbsOut$dW[-(1:burn),1], gibbsOut$dW[-(1:burn),2],
     xlab=expression(W[11]), ylab=expression(W[22]), pch='.')



###   Section 4.5.2: SUTSE models

inv <- read.table("Datasets/invest2.dat",col.names=c("Denmark","Spain"))
y <- ts(inv, frequency = 1, start = 1960) 

## prior hyperparameters
delta0 <- delta2 <- 3; delta1 <- 100
V0 <- (delta0-2) *diag(c(10^2, 500^2))
Wmu0 <- (delta1-2) * diag(0.01^2,2)
Wbeta0 <- (delta2 -2) * diag(c(5^2, 100^2))

## Gibbs sampling
MC <- 30000
TT <- nrow(y) 
gibbsTheta <- array(0, dim=c(TT+1,4, MC-1)) 
gibbsV <- array(0, dim=c(2,2, MC))
gibbsWmu <- array(0, dim=c(2,2, MC)) 
gibbsWbeta <- array(0, dim=c(2,2, MC)) 
mod <- dlm(FF = matrix(c(1,0),nrow=1) %x% diag(2),
           V = diag(2),
           GG = matrix(c(1,0,1,1),2,2) %x% diag(2),
           W = bdiag(diag(2), diag(2)),
           m0 = c(inv[1,1], inv[1,2],0,0),
           C0 = diag(x = 1e7, nrow = 4))
# starting values 
mod$V <- gibbsV[,,1] <- V0/(delta0-2)
gibbsWmu[,,1] <- Wmu0/(delta1-2)
gibbsWbeta[,,1] <- Wbeta0/(delta2-2)
mod$W <- bdiag(gibbsWmu[,,1], gibbsWbeta[,,1])

set.seed(3420)
for(it in 1: (MC-1))
	{
   	# generate states - FFBS
	modFilt <- dlmFilter(y, mod, simplify=TRUE)
	gibbsTheta[,,it] <- theta <- dlmBSample(modFilt)
	# update V
	S <- crossprod(y-theta[-1,] %*% t(mod$FF)) + V0
	gibbsV[,,it+1] <- solve(rwishart(df=delta0+1+TT,p=2,Sigma=solve(S)))
	mod$V <- gibbsV[,,it+1]
	# update Wmu and Wbeta
   	theta.center <- theta[-1,]-(theta[-(TT+1),]  %*% t(mod$GG))
      SS1 <- crossprod(theta.center)[1:2,1:2] + Wmu0
   	SS2 <- crossprod(theta.center)[3:4,3:4] + Wbeta0
   	gibbsWmu[,,it+1] <- solve(rwishart(df=delta1+1+TT, Sigma=solve(SS1)))
   	gibbsWbeta[,,it+1] <- solve(rwishart(df=delta2+1+TT, Sigma=solve(SS2)))
	mod$W <- bdiag(gibbsWmu[,,it+1], gibbsWbeta[,,it+1])   	
	}

## MCMC diagnostics
burn<- 1:20000
#traces
#ts.plot(sqrt(gibbsV[1,1,-burn]), xlab="iteration", ylab=expression(sigma[1]))
#ts.plot(sqrt(gibbsV[2,2,-burn]), xlab="iteration", ylab=expression(sigma[2]))
#ts.plot(gibbsV[2,1,-burn], xlab="iteration", ylab=expression(sigma[1][2]))
# V
par(mar = c(2, 4, 1, 1) + 0.1, cex = 1)
par(mfrow=c(3,2))
plot(ergMean(sqrt(gibbsV[1,1, -burn])),type="l", main="",cex.lab=1.5, 
     	ylab=expression(sigma[1]), xlab="MCMC iteration")
acf(sqrt(gibbsV[1,1,-burn]),  main="")
plot(ergMean(sqrt(gibbsV[2,2, -burn])),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[2]), xlab="MCMC iteration")
acf(sqrt(gibbsV[2,2,-burn]),  main="")
plot(ergMean(gibbsV[2,1, -burn]),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[1][2]), xlab="MCMC iteration")
acf(gibbsV[2,1,-burn],  main="")
# Wmu
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.8)
par(mfrow=c(3,2))
plot(ergMean(sqrt(gibbsWmu[1,1, -burn])),type="l", main="", cex.lab=1.5, 
	ylab=expression(sigma[mu][1]), xlab="MCMC iteration")
acf(sqrt(gibbsWmu[1,1,-burn]),  main="")
plot(ergMean(sqrt(gibbsWmu[2,2, -burn])),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[mu][2]), xlab="MCMC iteration")
acf(sqrt(gibbsWmu[2,2,-burn]),  main="")
plot(ergMean(gibbsWmu[2,1, -burn]),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[mu][1][2]), xlab="MCMC iteration")
acf(gibbsWmu[2,1,-burn],  main="")
# Wbeta
par(mar = c(2, 4, 1, 1) + 0.1, cex =1)
par(mfrow=c(3,2))
plot(ergMean(sqrt(gibbsWbeta[1,1, -burn])),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[beta][1]), xlab="MCMC iteration")
acf(sqrt(gibbsWbeta[1,1,-burn]),  main="")
plot(ergMean(sqrt(gibbsWbeta[2,2, -burn])),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[beta][2]), xlab="MCMC iteration")
acf(sqrt(gibbsWbeta[2,2,-burn]),  main="")
plot(ergMean(gibbsWbeta[2,1, -burn]),type="l", main="",cex.lab=1.5, 
	ylab=expression(sigma[beta][1][2]), xlab="MCMC iteration")
acf(gibbsWbeta[2,1,-burn],  main="")

## Synthesis of MCMC output
estV <- cbind(mcmcMean(gibbsV[1,1,-burn]),mcmcMean(gibbsV[2,2,-burn]),
		mcmcMean(gibbsV[2,1,-burn])); estV
estWmu <- cbind(mcmcMean(gibbsWmu[1,1,-burn]),
           	mcmcMean(gibbsWmu[2,2,-burn]),
           	mcmcMean(gibbsWmu[2,1,-burn])); round(estWmu,4)
estWbeta <- cbind(mcmcMean(gibbsWbeta[1,1,-burn]),
          	mcmcMean(gibbsWbeta[2,2,-burn]),
          	mcmcMean(gibbsWbeta[2,1,-burn])); estWbeta

## Plot: data and estimated level 
meanGibbsTheta <- apply(gibbsTheta[,,-burn],c(1,2), mean)
qinf <- function(x){quantile(x, prob=.025)}
qinfGibbsTheta <- apply(gibbsTheta[,,-burn],c(1,2), qinf)
qsup <- function(x){quantile(x, prob=.975)}
qsupGibbsTheta <- apply(gibbsTheta[,,-burn],c(1,2), qsup)
par(mar = c(2, 4, 1, 1) + 0.1, cex=1)
par(mfrow=c(2,1))
require(zoo)
plot(as.zoo(y[,1]), main = "", xlab = "year", cex.lab = 0.7,
    	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	col="darkgray", ylab="Investments - Denmark", type="o", pch=20)
lines(as.zoo(ts(meanGibbsTheta[-1,1], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="o", pch=20,lty=3,cex=.8)
lines(as.zoo(ts(qinfGibbsTheta[-1,1], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="l", lty=2)
lines(as.zoo(ts(qsupGibbsTheta[-1,1], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="l", lty=2)
plot(as.zoo(y[,2]), main = "", xlab = "year", cex.lab = 0.7,
     	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	col="darkgray", ylab="Investments - Spain", type="o", pch=20)
lines(as.zoo(ts(meanGibbsTheta[-1,2], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="o", pch=20,lty=3,cex=.8)
lines(as.zoo(ts(qinfGibbsTheta[-1,2], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="l", lty=2)
lines(as.zoo(ts(qsupGibbsTheta[-1,2], freq=frequency(y), start=start(y))), 
	oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1), 
	type="l", lty=2)



### Section 4.5.3 - A model for outliers and structural breaks 
###                 (lambda-omega_t model)

## require function "dlmGibbsDIGt"

par(mar=c(5.1,4.1,1,1))
curve(qnorm(0.9, sd=(1.642*x)^(-0.5)), 0.03, 1.1, xlab=expression(omega),
      ylab="90th percentile")
points(1,qnorm(0.9, sd=(1.642)^(-0.5)), pch=4, cex=3)

y <- log(UKgas)
set.seed(4521)
MCMC <- 10500
gibbsOut <- dlmGibbsDIGt(y, mod = dlmModPoly(2) + dlmModSeas(4),
                         A_y = 10000, B_y = 10000, p = 3, n.sample = MCMC,
                         save.states = TRUE, thin = 2)

## output analysis   
burn <- 1 : 500
nuRange <- c(1 : 10, seq(20, 100, by = 10))
omega_y <- ts(colMeans(gibbsOut$omega_y[-burn, ]),
              start = start(y), freq=4)
omega_theta <- ts(apply(gibbsOut$omega_theta[,, -burn], 1 : 2,
                        mean), start = start(y), freq = 4)
layout(matrix(c(1, 2, 3, 4), 4, 1, TRUE))

par(mar = c(3, 5, 1, 1) + 0.1)
plot(omega_y, type = "p", ylim = c(0, 1.2), pch = 16,
     xlab = "", ylab = expression(omega[list(y, t)]), cex.lab = 1.6)
abline(h = 1, lty = "dashed")
for (i in 1 : 3)
{ 
    plot(omega_theta[,i], ylim=c(0,1.2), pch = 16,
         type = "p", xlab = "",
         ylab = bquote(omega[list(theta, t * .(i))]), cex.lab = 1.6)
    abline(h = 1, lty = "dashed")
}

## identify the mild observational outlier
outlier <- which.min(omega_y)
time(y)[outlier]
omega_y[outlier]

## identify time of most extreme break in seasonal component 
a <- which.min(omega_theta[, 3])
time(y)[a]
omega_theta[a, 3]

## posterior inference for states
thetaMean <- ts(apply(gibbsTheta, 1 : 2, mean),
                start = start(y),
                freq = frequency(y))
LprobLim <- ts(apply(gibbsTheta, 1 : 2, quantile,
                     probs = 0.025),
               start = start(y), freq = frequency(y))
UprobLim <- ts(apply(gibbsTheta, 1 : 2, quantile,
                     probs = 0.975),
               start = start(y), freq = frequency(y))

par(mfrow = c(2, 1), mar = c(5.1, 4.1, 2.1, 2.1))
plot(thetaMean[, 1], xlab = "", ylab = "Trend")
lines(LprobLim[, 1], lty = 2); lines(UprobLim[, 1], lty = 2)
plot(thetaMean[, 3], xlab = "", ylab = "Seasonal", type = "o")
lines(LprobLim[, 3], lty = 2); lines(UprobLim[, 3], lty = 2)


###
### Section 4.6 -- Examples
###


### 4.6.1  Example: Estimating the output gap: Bayesian inference

## require function "gdpGibbs"

gdp <- log(read.table("Datasets/gdp5004.dat"))
gdp <- ts(gdp, start = 1950, frequency = 4)

set.seed(4521)
MCMC <- 6000
system.time(
            gibbsOut <- gdpGibbs(gdp, a=1, b=1000, n=MCMC, thin = 9)
            )

str(gibbsOut)
burn <- 1000

### Graphical convergence diagnostics
plot(gibbsOut$phi[-(1:burn),1], cex=0.5, xlab="", ylab="")
plot(gibbsOut$phi[-(1:burn),2], cex=0.5, xlab="", ylab="")
plot(gibbsOut$vars[-(1:burn),1], cex=0.5, xlab="", ylab="")
plot(gibbsOut$vars[-(1:burn),2], cex=0.5, xlab="", ylab="")
plot(gibbsOut$vars[-(1:burn),3], cex=0.5, xlab="", ylab="")
use <- MCMC - burn
from <- 0.05 * use
plot(ergMean(gibbsOut$phi[-(1:burn),1], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$phi[-(1:burn),2], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$vars[-(1:burn),1], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$vars[-(1:burn),2], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
plot(ergMean(gibbsOut$vars[-(1:burn),3], from), type="l",
     xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))
acf(gibbsOut$phi[-(1:burn),1], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$phi[-(1:burn),2], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$vars[-(1:burn),1], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$vars[-(1:burn),2], ci=0, ylim=c(0,1), ylab="", main="")
acf(gibbsOut$vars[-(1:burn),3], ci=0, ylim=c(0,1), ylab="", main="")

### Bayes estimates of unknown parameters
mcmcMean(gibbsOut$phi[-(1:burn),])
mcmcMean(gibbsOut$vars[-(1:burn),])


###  4.6.2 Example: Dynamic regression. 

## Data
yields <- read.table("Datasets/yields.dat")
y <- yields[1:192,3:19]  # leave out very short times to maturity
y <- as.matrix(y)
x <- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
p <- 3; m <- ncol(y) ; TT=nrow(y)
persp(x=x, z=t(y), theta=40, phi=30, expand=.5, col="lightgrey",
      ltheta=100,shade=0.75, 
      xlab="maturity (months)", ylab="time", zlab="yield")

## Cross-sectional model (Nelson and Siegel)
lambda <- 0.0609
h2 <- function(x){(1-exp(-lambda*x))/(lambda*x)}
h3 <- function(x){((1-exp(-lambda*x))/(lambda*x)) - exp(-lambda*x)}
X=cbind(rep(1,m), h2(x), h3(x))
# OLS estimates
betahat <- solve(crossprod(X), crossprod(X, t(y)))
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
plot.ts(t(betahat), lty=1:3, ylab="OLS estimates")
# Variance estimates
yfit <- t(X %*% betahat)
res <- (y-yfit)^2
s2 <- rowSums(res) / (m-p)
ts.plot(sqrt(s2))
# Dynamics of the OLS estimates
acf(t(betahat)) 
acf(diff(t(betahat)))

# Estimated regression curves (yields curves)
nelsonSiegel <- function(x, beta)
			{beta[1]+ beta[2]*h2(x)+ beta[3]*h3(x)}
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(2,2))
month=51
plot(x,y[month,], xlab="maturity(months)", ylab="yield (percent)", 
		ylim=c(8.9,9.8), main="3/31/89", cex=.6)
lines(x, nelsonSiegel(x, betahat[,month]))
month=55
plot(x,y[month,], xlab="maturity(months)", ylab="yield (percent)", 
	main="7/31/89", cex=.6)
lines(x, nelsonSiegel(x, betahat[,month]))
month=149
plot(x,y[month,], xlab="maturity(months)", ylab="yield (percent)", 
	main="5/30/97", cex=.6)
lines(x, nelsonSiegel(x, betahat[,month]))
month=164
plot(x,y[month,], xlab="maturity(months)", ylab="yield (percent)", 
	main="8/31/98", cex=.6)
lines(x, nelsonSiegel(x, betahat[,month]))

## Dynamic regression.  

mod <- dlm(m0=rep(0,p), C0=100 * diag(p), 
		FF=X, V=diag(m), GG=diag(p), W=diag(p))
##  Prior hyperparameters
#   ** Prior on the AR(1) coefficients: 
#      independent Gaussians, not resricted on the stationariety region
#             psi_j \sim N(0, 1)  , j=1,..,p(=3). 
##  ** Prior on the unknown variances
#      * Observation errors:
#            V_i \sim inv-gamma(a,b) , i=1,..., m(=17) 
#        taking a=3, b=.01, we have E(V)=0.005, sd(V)=.005
#      * State evolution errors
#            W_j \sim inv-gamma(aa,bb), j=1,..,p(=3); 
#        taking aa=3, bb=1, we have E(W)=0.5,  sd(W)=.5
psi0 <- 0; tau0 <- 1
shapeY <- 3; rateY <- .01
shapeTheta <- 3; rateTheta <- 1

## MCMC 
MC <- 10000  
gibbsTheta <- array(NA, dim=c(MC,TT+1,p))
gibbsPsi <- matrix(NA, nrow=MC, ncol=p)
gibbsV <- matrix(NA, nrow=MC, ncol=m)
gibbsW <- matrix(NA, nrow=MC, ncol=p)

# (0) Starting values -- as specified by mod 
set.seed(3420)
phi.init <- rnorm(3,psi0, tau0)
V.init <- 1/rgamma(1,shapeY,rateY)
W.init <- 1/rgamma(1,shapeTheta,rateTheta)
mod$GG=diag(phi.init)
mod$V=diag(rep(V.init, m))
mod$W=diag(rep(W.init, p))

## Gibbs sampling
for (i in 1:MC)
	{
	# generate the states by FFBS
	modFilt <- dlmFilter(y, mod, simplify=TRUE)
	theta <- dlmBSample(modFilt)
	gibbsTheta[i,,] <- theta
	# generate the W_j
	theta.center <- theta[-1,] - t(mod$GG %*% t(theta[-(TT+1),]))
	SStheta=apply((theta.center)^2, 2,sum)
	phiTheta=rgamma(p, shape=shapeTheta+TT/2, rate=rateTheta+SStheta/2)
	gibbsW[i,] <- 1/phiTheta
	mod$W <- diag(gibbsW[i,])
	# generate the V_i
	y.center <- y - t(mod$FF %*% t(theta[-1,]))
	SSy=apply((y.center)^2, 2,sum)
	gibbsV[i,] <- 1/rgamma(m, shape=shapeY+TT/2, rate=rateY+SSy/2)
	mod$V <- diag(gibbsV[i,])
	# generate the AR parameters psi_1, psi_2, psi_3
	psi.AR=rep(NA,3)
	for (j in 1:p)
		{tau= 1/((1/tau0) + phiTheta[j] * crossprod(theta[-(TT+1),j]))
		 psi=tau * (phiTheta[j] * t(theta[-(TT+1),j])  %*% theta[-1,j] + psi0/tau0)
		 psi.AR[j]=rnorm(1, psi, sd=tau^.5)
		}
	gibbsPsi[i,] <- psi.AR
	mod$GG <-diag(psi.AR)
	}

## MCMC output
burn <- 1000

## convergence diagnostics
# AR parameters
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,3))
for (j in 1:p)
	{
	 ts.plot(gibbsPsi[-(1:burn), j], xlab="iteration", ylab=paste("psi", as.character(j)))
	 ts.plot(cumsum((gibbsPsi[-(1:burn), j]))/(1:(MC-burn)), xlab="iteration",
    		ylab="ergMean")
	 acf(gibbsPsi[-(1:burn), j], main="")
	}
# Variances W_j
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,3))
for (j in 1:p)
	{
	 ts.plot(gibbsW[-(1:burn), j], xlab="iteration",ylab=paste("W", as.character(j)))
	 ts.plot(cumsum((gibbsW[-(1:burn), j]))/(1:(MC-burn)), xlab="iteration",
    		ylab="ergMean")
	 acf(gibbsW[-(1:burn), j], main="")
	}
# Variances V_i (only for a few i, i=1,5,17 say) 
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,3))
j=1	 
ts.plot(gibbsV[-(1:burn), j], xlab="iteration", ylab="V1 (3 months)")
ts.plot(cumsum((gibbsV[-(1:burn), j]))/(1:(MC-burn)), 
		xlab="iteration", ylab="ergMean")
acf(gibbsV[-(1:burn), j], main="")
j=5
ts.plot(gibbsV[-(1:burn), j], xlab="iteration", ylab="V5 (15 months)")
ts.plot(cumsum((gibbsV[-(1:burn), j]))/(1:(MC-burn)), 
		xlab="iteration", ylab="ergMean")
acf(gibbsV[-(1:burn), j], main="")
j=17
ts.plot(gibbsV[-(1:burn), j], xlab="iteration", ylab="V17(120 months)")
ts.plot(cumsum((gibbsV[-(1:burn), j]))/(1:(MC-burn)), 
		xlab="iteration", ylab="ergMean")
acf(gibbsV[-(1:burn), j], main="")

## Synthesis of MCMC output
round(mcmcMean(gibbsPsi[-(1:burn),]),4)
round(mcmcMean(gibbsW[-(1:burn),])^.5,4)
round(mcmcMean(gibbsV[-(1:burn),])^.5,4)
 
meanGibbsTheta <- apply(gibbsTheta[-burn,,], c(2,3), mean)
meanTheta=meanGibbsTheta[-1,]
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
plot.ts(meanTheta[,1], ylab="states smoothing and OLS estimates", xlab="time", ylim=c(-5,12),lwd=2)
lines(meanTheta[,2], col="darkgray",lwd=2)
lines(meanTheta[,3], col="gray",lwd=2)
lines(betahat[1,], lty=2)
lines(betahat[2,], lty=2)
lines(betahat[3,], lty=2)



###  4.6.3 Example: Factor models -- common stochastic trend.        

# Data
interestRate <- read.table ("Datasets/interestRates.dat", 
				col.names=c("Long","Short"))
y <- log(1+interestRate/100)
y <- ts(y, frequency = 52, start = 1971)
par(mar=c(2,4,1,1)+ 0.1, cex=0.7)
ts.plot(y, lty=c(1,2), col=c(1, "darkgray"))
legend ("topright",legend = c("mortgage rate (long rate)", 
	"federal funds rate (short rate)"), 
	col=c(1, "darkgray"), lty=c(1,2), bty="n")

# Prior hyperparameters
alpha0 <- 0; tau2 <- 16
expSigmaMu<- 0.01; varSigmaMu<- 1
a <- (expSigmaMu^2/varSigmaMu)+2; b <- expSigmaMu*(a-1)
delta <- 3; 
V0=matrix(c(1, 0.5, 0.5, 4), byrow=T, nrow=2)  

# Gibbs sampling
MC <- 10000
n <- nrow(y) 
gibbsTheta <- array(0, dim=c(n+1,2,MC-1)) #MC-1 matrices (n+1)x2
gibbsV <-array(0, dim=c(2,2,MC)) # MC matrices 2x2
gibbsAlpha <- rep(0,MC) 
gibbsW <- rep(0,MC)
# model and starting values for the MCMC
mod <- dlmModPoly(2, dW=c(1,0), C0=100*diag(2)) 
gibbsAlpha[1] <- 0     
mod$FF <- rbind(c(1,0), c(gibbsAlpha[1],1))
mod$W[1,1] <- gibbsW[1] <- 1/rgamma(1, a, rate=b)    
mod$V <- gibbsV[,,1] <- V0 /(delta-2)
mod$GG <- diag(2)
# MCMC loop
set.seed(1562)
for(it in 1:(MC-1))
	{
	# generate state- FFBS
	modFilt <- dlmFilter(y, mod, simplify=TRUE)
   	gibbsTheta[,,it] <- theta <- dlmBSample(modFilt)
	# update alpha
	rho <- gibbsV[,,it][1,2]/(gibbsV[,,it][1,1]*gibbsV[,,it][2,2])^.5	
	tauT <- (gibbsV[,,it][2,2]*(1-rho^2)*tau2) /
	      (tau2 * sum(theta[-1,1]^2)+ (1-rho^2)*gibbsV[,,it][2,2])
      alphaT <- tauT * ((tau2/(gibbsV[,,it][2,2])^.5) * 
		sum(((y[,2]-theta[-1,2])/gibbsV[,,it][2,2]^.5 -
            rho* (y[,1]-theta[-1,1])/gibbsV[,,it][1,1]^.5) * 
            theta[-1,1])+ alpha0*(1-rho^2))/(tau2*(1-rho^2))
	mod$FF[2,1] <- gibbsAlpha[it+1] <- rnorm(1, alphaT, tauT^.5)
	# update sigma_mu
 	SSmu <- sum( diff(theta[,1])^2)
      mod$W[1,1] <- gibbsW[it+1] <- 1/rgamma(1, a+n/2, rate=b+SSmu/2)
	# update V
	S <- V0 + crossprod(y- theta[-1,] %*% t(mod$FF))
 	mod$V <- gibbsV[,,it+1] <- solve(rwishart(df=delta+1+ n, 
			Sigma=solve(S)))
	}

## MCMC diagnostics
# traces
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,2))
ts.plot(gibbsAlpha[-burn], xlab="iteration", ylab=expression(alpha))
ts.plot(gibbsW[-burn], xlab="iteration", ylab=expression(sigma[mu]))
ts.plot(sqrt(gibbsV[1,1,-burn]), xlab="iteration", ylab=expression(sigma[1]))
ts.plot(gibbsV[2,1,-burn], xlab="iteration", ylab=expression(sigma[1][2]))
ts.plot(sqrt(gibbsV[2,2,-burn]), xlab="iteration", ylab=expression(sigma[2]))
# alpha 
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,1))
hist(gibbsAlpha[-burn], prob=T, main="")
lines(density(gibbsAlpha[-burn]))
ts.plot(ergMean(gibbsAlpha[-burn]), xlab="iteration",ylab="ergMean")
acf(gibbsAlpha[-burn], main="")
# sigma_mu
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,1))
hist(gibbsW[-burn], prob=T, main="")
lines(density(gibbsW[-burn]))
ts.plot(ergMean(gibbsW[-burn]), xlab="iteration",ylab="ergMean")
acf(gibbsW[-burn], main="")
# V
par(mar = c(2, 4, 1, 1) + 0.1, cex = 0.6)
par(mfrow=c(3,2))
plot(ergMean(sqrt(gibbsV[1,1, -burn])),type="l", main="",
	ylab=expression(sigma[1]), xlab="MCMC iteration")
acf(sqrt(gibbsV[1,1,-burn]),  main="")
plot(ergMean(sqrt(gibbsV[2,2, -burn])),type="l", main="",
	ylab=expression(sigma[2]), xlab="MCMC iteration")
acf(sqrt(gibbsV[2,2,-burn]),  main="")
plot(ergMean(gibbsV[2,1, -burn]),type="l", main="",
	ylab=expression(sigma[1][2]), xlab="MCMC iteration")
acf(gibbsV[2,1,-burn],  main="")

## Synthesis of MCMC output
# alpha
round(mcmcMean(gibbsAlpha[-burn]),5)
round(quantile(gibbsAlpha[-burn], probs=c(0.05,0.95)),5)
# State variance 
round(mcmcMean(sqrt(gibbsW[-burn])),5)
round(quantile(sqrt(gibbsW[-burn]), probs=c(0.05,0.95)),5)
# Observation covariance matrix V
varV <- apply(gibbsV[,,-burn],c(1,2),mean); varV
round(mcmcMean(sqrt(gibbsV[1,1, -burn])),5) # estimate of sigma_1
round(quantile(sqrt(gibbsV[1,1, -burn]), probs=c(0.05,0.95)),5)
round(mcmcMean(sqrt(gibbsV[2,2, -burn])),5)  #estimate of  sigma_v2
round(quantile(sqrt(gibbsV[2,2, -burn]), probs=c(0.05,0.95)),5)
round(mcmcMean(gibbsV[2,1, -burn]),5)  #estimate of cov V12
round(quantile(gibbsV[2,1, -burn], probs=c(0.05,0.95)),5)

# Plot: states -- common trend 
meanGibbsTheta <- apply(gibbsTheta[,,-burn],c(1,2), mean)
par(mar=c(2,4,1,1)+ 0.1, cex=1)
ts.plot(y, col=c(1,"darkgray"), xlab="year", ylab="", ylim=c(0,.2)) 
lines(ts(meanGibbsTheta[-1,1], freq=52, start=1971),lty=3, lwd=2)
legend ("topright",legend = c("long","short","Common trend"), lty=c(1,1,3),
		col=c(1,"darkgray",1), bty="n")
