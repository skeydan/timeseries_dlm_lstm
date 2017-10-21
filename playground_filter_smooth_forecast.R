library(dlm)
library(ggplot2)
library(forecast)

#y <- as.ts(c(1:10))
y <- AirPassengers
y <- window(AirPassengers, start = 1958)

model <- dlmModPoly(order = 1, dV = 6000, dW = 6000)
unlist(model)


################ filtering #################

f <- dlmFilter(y, model)
str(f,1)
# predictive distribution of the state
f$a
# filtered distribution of the state, starts 1 early
f$m
# one-step ahead forecasts
f$f

n <- length(y)
# variance on the last state
dlmSvd2var(f$U.C[[n+1]], f$D.C[n+1, ])


plot(y, type='o', col = c("darkgrey"), xlab = "", ylab = "Level")
mod1 <- dlmModPoly(order = 1, dV = 6000, dW = 6000)
filt1 <- dlmFilter(y, mod1)
lines(dropFirst(filt1$m), lty = "longdash", col="red")
mod2 <- dlmModPoly(order = 1, dV = 12000, dW = 2000)
filt2 <- dlmFilter(y, mod2)
lines(dropFirst(filt2$m), lty = "dotdash", col="green")
mod3 <- dlmModPoly(order = 1, dV = 2000, dW = 12000)
filt3 <- dlmFilter(y, mod3)
lines(dropFirst(filt3$m), lty = "dotted", col="blue")


################ smoothing #################

s <- dlmSmooth(f)
str(s, 1)
drop(dlmSvd2var(s$U.S[[n + 1]], s$D.S[n + 1,]))
drop(dlmSvd2var(s$U.S[[n / 2 + 1]], s$D.S[n / 2 + 1,]))

hwid <- qnorm(0.025, lower = FALSE) * sqrt(unlist(dlmSvd2var(s$U.S, s$D.S)))
smooth <- cbind(s$s, as.vector(s$s) + hwid %o% c(-1, 1))
plot(dropFirst(smooth), plot.type = "s", type = "l", lty = c(1, 5, 5), ylab = "Level", xlab = "",
     ylim = range(y))
lines(y, type = "o", col = "darkgrey")

########################################

model <- dlm(
  m0 = rep(0, 4),
  C0 = 1e8 * diag(4),
  FF = matrix(c(1, 1, 0, 0), nr = 1),
  V = 1e-3,
  GG = bdiag(matrix(1),
             matrix(
               c(-1,-1,-1, 1, 0, 0, 0, 1, 0),
               nr = 3, byrow = TRUE
             )),
  W = diag(c(771.35, 86.48, 0, 0), nr = 4)
)
plot(y, xlab = "", ylab = "", type = 'o', col = "darkgrey")
### Filter
filt <- dlmFilter(y, model)
lines(dropFirst(filt$m[, 1]), lty = "dotdash")
### Smooth
smooth <- dlmSmooth(filt)
lines(dropFirst(smooth$s[,1]), lty = "longdash")
### Seasonal component
plot(dropFirst(smooth$s[, 3]), type = 'o', xlab = "")
abline(h = 0)


################ forecasting #################

# one-step-ahead forecasts for different signal-to-noise ratios
a <- window(cbind(y, filt1$f, filt2$f))
plot(a[, 1], type = 'o', col = "darkgrey", xlab = "", ylab = "Level")
lines(a[, 2], lty = "longdash")
lines(a[, 3], lty = "dotdash")

# signal-to-noise ratio varies
# Nile dam example

mod0 <- dlmModPoly(order = 1, dV = 15100, dW = 1468)
X <- ts(matrix(mod0$W, nc = 1, nr = length(Nile)), start = start(Nile))
window(X, 1898, 1899) <- 12 * mod0$W
modDam <- mod0
modDam$X <- X
modDam$JW <- matrix(1, 1, 1)
damFilt <- dlmFilter(Nile, modDam)
mod0Filt <- dlmFilter(Nile, mod0)
a <- window(cbind(Nile, mod0Filt$f, damFilt$f),start = 1880, end = 1920)
plot(a[, 1], type = 'o', col = "darkgrey", xlab = "", ylab = "Level")
lines(a[, 2], lty = "longdash")
lines(a[, 3], lty = "dotdash")
abline(v = 1898, lty = 2)
qqnorm(residuals(damFilt, sd = FALSE))
qqline(residuals(damFilt, sd = FALSE))
tsdiag(damFilt)


# n-step ahead forecasts
model <- dlmModPoly(order = 1, dV = 6000, dW = 6000)
f <- dlmFilter(y, model)
fc <- dlmForecast(f, nAhead = 12, sampleNew = 6)
plot(y, xlim=c(1958, 1962))
names(fc)
invisible(lapply(fc$newObs, function(x) lines(x, col = "darkgrey", type = 'o', pch = 4)))
lines(fc$f, type = 'o', lwd = 2, pch = 16)
abline(v = mean(c(time(fc$f)[1], time(y)[length(y)])), lty = "dashed")


# signal-to-noise ratio and forecasting
# needs no data!
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