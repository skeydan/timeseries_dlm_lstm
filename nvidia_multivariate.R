library(ggplot2)
library(forecast)
library(dlm)
library(zoo)
library("quantmod")

amd <- getSymbols("AMD", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$AMD.Adjusted
autoplot(amd)
nvda <- getSymbols("NVDA", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$NVDA.Adjusted
autoplot(nvda)
smh <- getSymbols("SMH", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$SMH.Adjusted
autoplot(smh)
mu <- getSymbols("MU", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$MU.Adjusted
autoplot(mu)


df <- data.frame(index(amd), amd, nvda, mu, smh)
(tss <- read.zoo(df))
autoplot(tss) + facet_grid(Series ~ ., scales = "free_y")

other <- mu


###   Section 4.5.2: SUTSE models
df <- data.frame( mu, nvda)
inv <- read.table("Datasets/invest2.dat",col.names=c("Denmark","Spain"))
y <- ts(df)

## prior hyperparameters
delta0 <- delta2 <- 3; delta1 <- 100
V0 <- (delta0-2) *diag(c(10^2, 500^2))
Wmu0 <- (delta1-2) * diag(0.01^2,2)
Wbeta0 <- (delta2 -2) * diag(c(5^2, 100^2))

## Gibbs sampling
MC <- 5000
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
burn<- 1:2000
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

