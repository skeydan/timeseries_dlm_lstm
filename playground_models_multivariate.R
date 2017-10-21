library(dlm)
library(ggplot2)
library(forecast)


###Example - Denmark and Spain Investments
## - integrated random walk SUTSE 
invest <- ts(matrix(scan("Datasets/invest2.dat"), nc = 2, byrow = TRUE),
             start = 1960, names = c("Denmark", "Spain"))

# system errors of level and slope may be correlated (W)
# evolution of levels and slopes assumed to be uncorrelated (G, F)
mod <- dlmModPoly(2)
str(mod)
mod$FF <- mod$FF %x% diag(2)
mod$FF
mod$GG <- mod$GG %x% diag(2)
mod$GG
mod$W
W1 <- matrix(0, 2, 2)     # level covariance is 0 (integrated random walk)
W2 <- diag(c(49, 437266)) # slope covariance
W2[1, 2] <- W2[2, 1] <- 155 # !!
mod$W <- bdiag(W1, W2)
mod$W
V <- diag(c(72, 14353)) # noise covariance
V[1, 2] <- V[2, 1] <- 1018 # !!
mod$V <- V
mod$V
mod$m0 <- rep(0, 4)
mod$C0 <- diag(4) * 1e7
mod <- as.dlm(mod)
investFilt <- dlmFilter(invest, mod)
par(mar=c(3.1,4.1,1.1,2.1), cex=0.5)
require(zoo)
plot(as.zoo(invest), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5, type='o')

sdev <- residuals(investFilt)$sd
lwr <- investFilt$f + qnorm(0.25) * sdev
upr <- investFilt$f - qnorm(0.25) * sdev

par(mar=c(2, 3, 1, 0) + 0.1, cex=0.7)
plot(invest[,1], type='o', ylim=c(32, 247), xlab="", ylab="")
lines(window(investFilt$f[,1], start = start(invest) + c(2,0)), type='o', lty=2, pch=4)
lines(window(lwr[,1], start = start(invest) + c(2,0)), col="darkgrey")
lines(window(upr[,1], start = start(invest) + c(2,0)), col="darkgrey")
legend("topleft", inset = 0.05,
       legend=c("Observed", "One-step-ahead forecast", "50% prediction interval"),
       pch=c(1,4,-1), lty=c(1,2,1), col=c(rep("black", 2), "darkgrey"), bty='n') 

par(mar=c(2, 3, 1, 0) + 0.1, cex=0.7)
plot(invest[,2], type='o', xlab="", ylab="")
lines(window(investFilt$f[,2], start = start(invest) + c(2,0)), type='o', lty=2, pch=4)
lines(window(lwr[,2], start = start(invest) + c(2,0)), col="darkgrey")
lines(window(upr[,2], start = start(invest) + c(2,0)), col="darkgrey")
legend("topleft", inset = 0.05,
       legend=c("Observed", "One-step-ahead forecast", "50% prediction interval"),
       pch=c(1,4,-1), lty=c(1,2,1), col=c(rep("black", 2), "darkgrey"), bty='n') 


######
### Section 3.3.3: Seemingly unrelated regression models
######

### Dynamic CAPM - SUTSE
tmp <- ts(read.table("Datasets/P.dat",
                     header = TRUE),
          start = c(1978, 1), frequency = 12) * 100
y <- tmp[, 1 : 4] - tmp[, "RKFREE"]
colnames(y) <- colnames(tmp)[1 : 4]
market <- tmp[, "MARKET"] - tmp[, "RKFREE"]
rm("tmp")
m <- NCOL(y)
### Set up the model
CAPM <- dlmModReg(market)
CAPM
CAPM$FF <- CAPM$FF %x% diag(m)
CAPM$GG <- CAPM$GG %x% diag(m)
CAPM$JFF <- CAPM$JFF %x% diag(m)
CAPM$W <- CAPM$W %x% matrix(0, m, m) # intercepts are uncorrelated
CAPM$W[-(1 : m), -(1 : m)] <-       # correlated slope variance
  c(8.153e-07,  -3.172e-05, -4.267e-05, -6.649e-05,
    -3.172e-05, 0.001377,   0.001852,   0.002884,  
    -4.267e-05, 0.001852,   0.002498,   0.003884,  
    -6.649e-05, 0.002884,   0.003884,   0.006057)  
CAPM$V <- CAPM$V %x% matrix(0, m, m) # correlated noise variance
CAPM$V[] <- c(41.06,     0.01571, -0.9504, -2.328,
              0.01571, 24.23,     5.783,   3.376, 
              -0.9504,   5.783,   39.2,     8.145, 
              -2.328,    3.376,    8.145,  39.29) 
CAPM$m0 <- rep(0, 2 * m)
CAPM$C0 <- diag(1e7, nr = 2 * m)
CAPM

## Smooth
CAPMsmooth <- dlmSmooth(y, CAPM)

## plots
par(mar = c(3, 4, 1, 2) + 0.1, cex = 0.7)
plot(dropFirst(CAPMsmooth$s[, m + 1 : m]),
     lty = c("13", "6413", "431313", "B4"),
     plot.type = "s", xlab = "", ylab = "Beta")
abline(h = 1, col = "darkgrey")
legend("bottomright", legend = colnames(y), bty = "n",
       lty = c("13", "6413", "431313", "B4"), inset = 0.05)
