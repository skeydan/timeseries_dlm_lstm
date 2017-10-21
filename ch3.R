### R code for Chapter 3
require(dlm)

######
### Section 3.2.1: Trend models
######

### Signal-to-noise ratio in local level model
par(mfrow = c(2, 2), mar = c(2, 2, 3, 1) + 0.1, cex = 0.5)
r <- 10^(seq(-1,2,length = 4))
for (i in 1 : 4)
{
    mod <- dlmModPoly(1, dV = 1 / r[i], dW = 1, m0 = 5, C0 = 0)
    set.seed(2)
    sim <- dlmForecast(mod, nAhead = 50, sampleNew = 1)
    if (i == 1) rng <- range(sim$newObs[[1]])
    plot(sim$newStates[[1]], ylim = rng, type = "l",
         xlab = "", ylab = "", col = "darkgrey",
         main = bquote(log[1*0] ~~ r   ==  .(r[i])))
    lines(sim$newObs[[1]], type = "o", cex = 0.5)
}

### local level model
lakeSup <- ts(read.table("Datasets/lakeSuperior.dat", skip = 3,
                         colClasses = "numeric")[, 2], start = 1900)
par(mfrow = c(1, 1))
par(mar=c(3.1,4.1,1.1,2.1), cex=0.7)
plot(lakeSup, type='o', xlab="", ylab="Level")
modLSup <- dlmModPoly(1, dV = 9.465, dW = 0.121)
lSupFilt <- dlmFilter(lakeSup, modLSup)
lSupSmooth <- dlmSmooth(lSupFilt)
### filtering & smoothing estimates
par(mfrow = c(2, 2), mar = c(3, 4, 1, 1) + 0.1, cex = 0.5)
plot(lakeSup, xlab = "", ylab = "(a)", col = "darkgrey")
lines(dropFirst(lSupFilt$m), lty = "longdash")
sdFilt <- with(lSupFilt, ts(sqrt(unlist(dlmSvd2var(U.C, D.C))),
               start = start(lakeSup)))
lns <- cbind(dropFirst(lSupFilt$m), 
             qnorm(0.95, sd = sdFilt[-1]) %o% c(-1, 1) +
             lSupFilt$m[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash")
plot(dropFirst(sdFilt), xlab = "", ylab = "(b)")
plot(lakeSup, xlab = "", ylab = "(c)", col = "darkgrey")
lines(dropFirst(lSupSmooth$s), lty = "longdash")
sdSmooth <- with(lSupSmooth, ts(sqrt(unlist(dlmSvd2var(U.S, D.S))),
               start = start(lakeSup)))
lns <- cbind(dropFirst(lSupSmooth$s), 
             qnorm(0.95, sd = sdSmooth[-1]) %o% c(-1, 1) +
             lSupSmooth$s[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash")
plot(dropFirst(sdSmooth), xlab = "", ylab = "(d)")

### Diagnostic plots
res <- residuals(lSupFilt, sd=FALSE)
plot(res,type='h'); abline(h=0)
acf(res)
qqnorm(res); qqline(res)
### Formal diagnostic tests
shapiro.test(res)
Box.test(res, lag=20, type="Ljung")
sapply(1 : 20, function(i)
       Box.test(res, lag = i, type = "Ljung-Box")$p.value)

### Holt-Winters (local level)
par(mfrow = c(1, 1))
par(mar=c(3.1,4.1,1.1,2.1), cex=0.7)
HWout <- HoltWinters(lakeSup, gamma = FALSE, beta = FALSE)
plot(dropFirst(lSupFilt$f), lty = "dashed",
     xlab = "", ylab = "")
lines(HWout$fitted[, "level"])
leg <- c("Holt-Winters", "Local level DLM")
legend("topleft", legend = leg, bty = "n",
       lty = c("solid", "dashed"))


### linear growth model
invSpain <- ts(read.table("Datasets/invest2.dat", colClasses = "numeric")[,2],
               start = 1960)
mod1 <- dlmModPoly(dV = 10, dW = c(102236, 321803))
mod1Filt <- dlmFilter(invSpain, mod1)
fut1 <- dlmForecast(mod1Filt, n=5)
mod2 <- dlmModPoly(dV = 10, dW = c(0, 515939))
mod2Filt <- dlmFilter(invSpain, mod2)
fut2 <- dlmForecast(mod2Filt, n=5)

mean(abs(mod1Filt$f - invSpain))
mean(abs(mod2Filt$f - invSpain))
mean((mod1Filt$f - invSpain)^2)
mean((mod2Filt$f - invSpain)^2)
mean(abs(mod1Filt$f - invSpain) / invSpain)
mean(abs(mod2Filt$f - invSpain) / invSpain)

sqrt(sum((mod1Filt$f - invSpain)[-(1:5)]^2) /
     sum(diff(invSpain[-(1:4)])^2))
sqrt(sum((mod2Filt$f - invSpain)[-(1:5)]^2) /
     sum(diff(invSpain[-(1:4)])^2))

par(mar=c(2,3,1,0) + 0.1, cex=0.7)
plot(invSpain, xlim=c(1960, 2005), ylim=c(2500, 22000),
     xlab="", ylab="Investments", type='o', col = "darkgrey")
lines(dropFirst(mod1Filt$f), type='o', lty="431313", pch=3)
lines(dropFirst(mod2Filt$f), type='o', lty="22848222", pch=4)
lines(fut1$f, lty="dashed")
lines(fut2$f, lty="dotted")
legend("topleft", bty = "n",
       legend=c("observed", "one-step-ahead forecast - model 1",
       "one-step-ahead forecast - model 2", "5 years forecast - model 1",
       "5 years forecast - model 2"),
       lty=c("solid", "431313", "22848222", "dashed", "dotted"), pch=c(1,3,4,-1,-1),
       col=c("darkgrey", rep("black", 4)), inset = 0.05)


### Seasonal factor models
mod <- dlmModSeas(frequency = 4, dV = 3.5, dW = c(4.2, 0, 0))

######
### Section 3.2.3: Fourier form of seasonal models
#####

### Plot harmonic functions
n <- 12 # even
omega <- 2 * pi / n
par(mfrow = c(n - 1, 1), mar = c(0, 3.1, 0, 3.1),
    oma = c(3, 0, 2, 0), pch = 16)
for (i in 1:(n/2 - 1)) {
    curve(cos(x * i * omega), 0, n, ylim = c(-1.2, 1.2),
          xlim = c(-0.2, n + 0.2), ylab = "", axes = FALSE)
    points(1:n, cos(i * omega * 1:n), cex = 1.5)
    axis(2, at = c(-1, 0, 1)); abline(h = 0, col = "grey")
    if (i == 1) {
        usr <- par("usr")
        segments(x0 = c(rep(usr[1], 2), usr[2]),
                 y0 = c(usr[4] - 2 * (n - 1), rep(usr[4], 2)),
                 x1 = c(usr[1], rep(usr[2], 2)),
                 y1 = c(rep(usr[4], 2), usr[4] - 2 * (n - 1)),
                 xpd = NA)
    }
    curve(sin(x * i * omega), 0, n, ylim = c(-1.2, 1.2),
          xlim = c(-0.2, n + 0.2), ylab = "", axes = FALSE)
    points(1:n, sin(i * omega * 1:n), cex = 1.5)
    axis(4, at = c(-1, 0, 1)); abline(h = 0, col = "grey")
}
curve(cos(x * (n/2) * omega), 0, n, ylim = c(-1.2, 1.2),
      xlim = c(-0.2, n + 0.2), ylab = "", axes = FALSE)
points(1:n, rep(c(-1,1), n/2), cex = 1.5)
axis(1); axis(2, at = c(-1, 0, 1)); abline(h = 0, col = "grey")
usr <- par("usr")
segments(x0 = c(rep(usr[1], 2), usr[2]),
         y0 = c(usr[3] + 2 * (n - 1), rep(usr[3], 2)),
         x1 = c(usr[1], rep(usr[2], 2)),
         y1 = c(rep(usr[3], 2), usr[3] + 2 * (n - 1)),
         xpd = NA)

### monthly data: full versus reduced model
mod1 <- dlmModTrig(s = 12, dV = 5.1118, dW = 0) +
    dlmModPoly(1, dV = 0, dW = 81307e-3)
smoothTem1 <- dlmSmooth(nottem, mod1)
mod2 <- dlmModTrig(s = 12, q = 2, dV = 5.1420, dW = 0) +
    dlmModPoly(1, dV = 0, dW = 81942e-3)
smoothTem2 <- dlmSmooth(nottem, mod2)
## MAPE full model
mean(abs(residuals(dlmFilter(nottem, mod1),
                   type = "raw", sd = FALSE)) / nottem)
## MAPE reduced model
mean(abs(residuals(dlmFilter(nottem, mod2),
                   type = "raw", sd = FALSE)) / nottem)
## plot harmonics
plot(ts(smoothTem1$s[2 : 13, c(1, 3, 5, 7, 9, 11)],
        names = paste("S", 1 : 6, sep = "_")),
     oma.multi = c(2, 0, 1, 0), pch = 16, nc = 1,
     yax.flip = TRUE, type = 'o', xlab = "", main = "")

######
### Section 3.2.4: General periodic componens
######

### Sunspots data
mod <- dlmModTrig(q = 2, tau = 130.51, dV = 0,
                  dW = rep(c(1765e-2, 3102e-4), each = 2)) +
        dlmModPoly(1, dV = 0.7452, dW = 0.1606)

sspots <- sqrt(sunspots)
sspots.smooth <- dlmSmooth(sspots, mod)
y <- cbind(sspots,
           tcrossprod(dropFirst(sspots.smooth$s[, c(1, 3, 5)]),
                      matrix(c(0, 0, 1, 1, 1, 0), nr = 2,
                             byrow = TRUE)))
colnames(y) <- c("Sunspots", "Level", "Periodic")
plot(y, yax.flip = TRUE, oma.multi = c(2, 0, 1, 0))

######
### Section 3.2.6: Example: estimating the output gap
######
gdp <- read.table("Datasets/gdp5004.dat")
gdp <- ts(gdp, frequency = 4, start = 1950)
Lgdp <- log(gdp)
par(mfrow = c(1, 1))
par(mar=c(3.1,4.1,1.1,2.1), cex=0.5)
plot(Lgdp, xlab = "", ylab = "log US GDP", main="")
### MLE - no stationarity constraints
level0 <- Lgdp[1]
slope0 <- mean(diff(Lgdp))
buildGap <- function(u) {
    trend <- dlmModPoly(dV = 1e-7, dW = exp(u[1:2]),
                        m0 = c(level0, slope0), C0 = 2 * diag(2))
    gap <- dlmModARMA(ar = u[4:5], sigma2 = exp(u[3]))
    return(trend + gap)}
init <- c(-3, -1, -3, .4, .4)
outMLE <- dlmMLE(Lgdp, init, buildGap)
dlmGap <- buildGap(outMLE$par)
sqrt(diag(W(dlmGap))[1:3])
GG(dlmGap)[3:4, 3]
## check stationarity of MLE solution
Mod(polyroot(c(1, -GG(dlmGap)[3:4, 3])))
plot(ARMAacf(ar = GG(dlmGap)[3:4, 3], lag.max = 20),
     ylim = c(0, 1), ylab = "acf", type = "h")
### Smoothing estimates
gdpSmooth <- dlmSmooth(Lgdp, dlmGap)
par(mar=c(3.1,4.1,1.1,2.1), cex=0.5)
plot(cbind(Lgdp, dropFirst(gdpSmooth$s[, 1])), plot.type = "single",
     xlab = "", ylab = "Log GDP", lty = c("longdash", "solid"),
     col = c("darkgrey", "black"))
plot(dropFirst(gdpSmooth$s[, c(1, 3)]), oma.multi = c(2,0,1,0),
     mar.multi = c(0, 2.1, 0, 2.1), cex.axis = 0.5,
     ann = FALSE, yax.flip = TRUE)
plot(dropFirst(gdpSmooth$s[, c(1, 3)]), ann = FALSE, yax.flip = TRUE)
### MLE - imposing stationarity constraints
buildgapr <- function(u)
{
    trend <- dlmModPoly(dV = 0.000001,
                        dW = c(exp(u[1]), exp(u[2])),
                        m0 = c(Lgdp[1], mean(diff(Lgdp))),
                        C0 = 2*diag(2))
    gap <- dlmModARMA(ar = ARtransPars(u[4 : 5]),
                      sigma2 = exp(u[3]))
    return(trend + gap)
}
init <- c(-3, -1, -3, .4, .4)
outMLEr <- dlmMLE(Lgdp, init, buildgapr)
outMLEr$value
parMLEr <- c(exp(outMLEr$par[1 : 3])^.5,
             ARtransPars(outMLEr$par[4 : 5]))
round(parMLEr, 4)

Mod(polyroot(c(1,-ARtransPars(outMLEr$par[4:5]))))
plot(ARMAacf(ar=ARtransPars(outMLEr$par[4:5]),lag.max=10),
    ylim=c(0,1), ylab="acf")
modr <- buildgapr(outMLEr$par)
outFr <- dlmFilter(Lgdp,modr)
outSr <- dlmSmooth(outFr)
ts.plot(cbind(Lgdp,outSr$s[-1,1]),col=1:2)
plot.ts(outSr$s[-1,c(1,3)],main="",ann=F,yax.flip=TRUE)



######
## Section 3.2.7: regression models dynamic 
######

###Example - Capital Asset Pricing Model
capm <- read.table("Datasets/P.dat", header = TRUE)
capm.ts <- ts(capm, start = c(1978, 1), frequency = 12)
colnames(capm)
par(cex = 0.5)
require(zoo)
plot(as.zoo(capm.ts), main = "", xlab = "", cex.lab = 0.7,
     oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1))
IBM <- capm.ts[, "IBM"] - capm.ts[, "RKFREE"]
x <- capm.ts[, "MARKET"] - capm.ts[, "RKFREE"]
outLM <- lm(IBM ~ x)
outLM$coef
summary(outLM)
acf(outLM$res)
qqnorm(outLM$res)
### static CAPM
mod <- dlmModReg(x, dV = 0.00254, m0 = c(0, 1.5),
                 C0 = diag(c(1e+07, 1)))
outF <- dlmFilter(IBM, mod)
outF$m[1 + length(IBM), ]

### dynamic CAPM
buildCapm <- function(u) {
    dlmModReg(x, dV = exp(u[1]), dW = exp(u[2 : 3]))
}
outMLE <- dlmMLE(IBM, parm = rep(0, 3), buildCapm)
exp(outMLE$par)
outMLE$value
mod <- buildCapm(outMLE$par)
outS <- dlmSmooth(IBM, mod)
require(zoo)
plot(as.zoo(dropFirst(outS$s)), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5)

######
### Section 3.3.2: Seemingly unrelated time series equations
######

###Example - Denmark and Spain Investments
## - integrated random walk SUTSE 
invest <- ts(matrix(scan("Datasets/invest2.dat"), nc = 2, byrow = TRUE),
             start = 1960, names = c("Denmark", "Spain"))
mod <- dlmModPoly(2)
mod$FF <- mod$FF %x% diag(2)
mod$GG <- mod$GG %x% diag(2)
W1 <- matrix(0, 2, 2)
W2 <- diag(c(49, 437266))
W2[1, 2] <- W2[2, 1] <- 155
mod$W <- bdiag(W1, W2)
V <- diag(c(72, 14353))
V[1, 2] <- V[2, 1] <- 1018
mod$V <- V
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
CAPM$FF <- CAPM$FF %x% diag(m)
CAPM$GG <- CAPM$GG %x% diag(m)
CAPM$JFF <- CAPM$JFF %x% diag(m)
CAPM$W <- CAPM$W %x% matrix(0, m, m)
CAPM$W[-(1 : m), -(1 : m)] <-
    c(8.153e-07,  -3.172e-05, -4.267e-05, -6.649e-05,
      -3.172e-05, 0.001377,   0.001852,   0.002884,  
      -4.267e-05, 0.001852,   0.002498,   0.003884,  
      -6.649e-05, 0.002884,   0.003884,   0.006057)  
CAPM$V <- CAPM$V %x% matrix(0, m, m)
CAPM$V[] <- c(41.06,     0.01571, -0.9504, -2.328,
               0.01571, 24.23,     5.783,   3.376, 
              -0.9504,   5.783,   39.2,     8.145, 
              -2.328,    3.376,    8.145,  39.29) 
CAPM$m0 <- rep(0, 2 * m)
CAPM$C0 <- diag(1e7, nr = 2 * m)
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
