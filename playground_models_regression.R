library(dlm)
library(ggplot2)
library(forecast)

###Example - Capital Asset Pricing Model
capm <- read.table("Datasets/P.dat", header = TRUE)
capm.ts <- ts(capm, start = c(1978, 1), frequency = 12)
colnames(capm.ts)

par(cex = 0.5)
require(zoo)
plot(as.zoo(capm.ts), main = "", xlab = "", cex.lab = 0.7, oma = c(2, 0, 1, 0), mar = c(0, 4.1, 0, 1.1))

IBM <- capm.ts[, "IBM"] - capm.ts[, "RKFREE"]
x <- capm.ts[, "MARKET"] - capm.ts[, "RKFREE"]

# lm
outLM <- lm(IBM ~ x)
outLM$coef
summary(outLM)
acf(outLM$res)
qqnorm(outLM$res)

### static CAPM
# dV is square of residual standard error from lm
# m0 and C0 are priors for mean and variance of alpha & beta
mod <- dlmModReg(x, dV = 0.00254, m0 = c(0, 1.5), C0 = diag(c(1e+07, 1)))
mod$W
mod$FF
mod$GG
mod$JFF #!!! not null
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
mod$FF
mod$GG
mod$W
outS <- dlmSmooth(IBM, mod)
require(zoo)
# development of coefficients
plot(as.zoo(dropFirst(outS$s)), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5)

