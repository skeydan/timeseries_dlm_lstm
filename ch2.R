###
rw <- dlm(m0 = 0, C0 = 10, FF = 1, V = 1.4, GG = 1, W = 0.2)
unlist(rw)
lg <- dlm(FF = matrix(c(1, 0), nr = 1),
          V = 1.4,
          GG = matrix(c(1, 0, 1, 1), nr = 2),
          W = diag(c(0, 0.2)),
          m0 = rep(0, 2),
          C0 = 10 * diag(2))
lg
is.dlm(lg)
V(lg) <- 0.8
W(lg)[2, 2] <- 0.5
V(lg)
W(lg)

###
x <- rnorm(100) # covariates
dlr <- dlm(FF = matrix(c(1, 0), nr = 1),
           V = 1.3,
           GG = diag(2),
           W = diag(c(0.4, 0.2)),
           m0 = rep(0, 2), C0 = 10 * diag(2),
           JFF = matrix(c(0, 1), nr = 1),
           X = x)
dlr

###
JV(rw) <- 1
is.dlm(rw)
dlm(rw)
X(rw) <- rep(c(0.75, 1.25), c(10, 20))
rw <- dlm(rw)
V(rw)

###
NilePoly <- dlmModPoly(order = 1, dV = 15100, dW = 1468)
unlist(NilePoly)
NileFilt <- dlmFilter(Nile, NilePoly)
str(NileFilt, 1)
n <- length(Nile)
attach(NileFilt)
dlmSvd2var(U.C[[n + 1]], D.C[n + 1, ])

###
plot(Nile, type='o', col = c("darkgrey"),
     xlab = "", ylab = "Level")
mod1 <- dlmModPoly(order = 1, dV = 15100, dW = 755)
NileFilt1 <- dlmFilter(Nile, mod1)
lines(dropFirst(NileFilt1$m), lty = "longdash")
mod2 <- dlmModPoly(order = 1, dV = 15100, dW = 7550)
NileFilt2 <- dlmFilter(Nile, mod2)
lines(dropFirst(NileFilt2$m), lty = "dotdash")
leg <- c("data", paste("filtered,  W/V =",
                       format(c(W(mod1) / V(mod1),
                                W(mod2) / V(mod2)))))
legend("bottomright", legend = leg,
       col=c("darkgrey", "black", "black"),
       lty = c("solid", "longdash", "dotdash"),
       pch = c(1, NA, NA), bty = "n")

###
NileSmooth <- dlmSmooth(NileFilt)
str(NileSmooth, 1)
attach(NileSmooth)
drop(dlmSvd2var(U.S[[n + 1]], D.S[n + 1,]))
drop(dlmSvd2var(U.C[[n + 1]], D.C[n + 1,]))
drop(dlmSvd2var(U.S[[n / 2 + 1]], D.S[n / 2 + 1,]))
drop(dlmSvd2var(U.C[[n / 2 + 1]], D.C[n / 2 + 1,]))
hwid <- qnorm(0.025, lower = FALSE) *
    sqrt(unlist(dlmSvd2var(U.S, D.S)))
smooth <- cbind(s, as.vector(s) + hwid %o% c(-1, 1))
plot(dropFirst(smooth), plot.type = "s", type = "l",
     lty = c(1, 5, 5), ylab = "Level", xlab = "",
     ylim = range(Nile))
lines(Nile, type = "o", col = "darkgrey")
legend("bottomleft", col = c("darkgrey", rep("black", 2)),
       lty = c(1, 1, 5), pch = c(1, NA, NA), bty = "n",
       legend = c("data", "smoothed level",
       "95% probability limits"))

###
expd <- ts(read.table("Datasets/qconsum.dat", skip = 4,
                      colClasses = "numeric")[, 1],
           start = c(1957, 1), frequency = 4)
expd.dlm <- dlm(m0 = rep(0,4), C0 = 1e8 * diag(4),
                FF = matrix(c(1, 1, 0, 0), nr = 1),
                V = 1e-3,
                GG = bdiag(matrix(1),
                matrix(c(-1, -1, -1, 1, 0, 0, 0, 1, 0),
                       nr = 3, byrow = TRUE)),
                W = diag(c(771.35, 86.48, 0, 0), nr = 4))
plot(expd, xlab = "", ylab = "Expenditures", type = 'o',
     col = "darkgrey")
### Filter
expdFilt <- dlmFilter(expd, expd.dlm)
lines(dropFirst(expdFilt$m[, 1]), lty = "dotdash")
### Smooth
expdSmooth <- dlmSmooth(expdFilt)
lines(dropFirst(expdSmooth$s[,1]), lty = "longdash")
legend("bottomright", col = c("darkgrey", rep("black", 2)),
       lty = c("solid", "dotdash", "longdash"),
       pch = c(1, NA, NA), bty = "n",
       legend = c("data", "filtered level", "smoothed level"))
### Seasonal component
plot(dropFirst(expdSmooth$s[, 3]), type = 'o', xlab = "",
     ylab = "Expenditure - Seasonal component")
abline(h = 0)

###
a <- window(cbind(Nile, NileFilt1$f, NileFilt2$f),
            start = 1880, end = 1920)
plot(a[, 1], type = 'o', col = "darkgrey",
     xlab = "", ylab = "Level")
lines(a[, 2], lty = "longdash")
lines(a[, 3], lty = "dotdash")
leg <- c("data", paste("one-step-ahead forecast,  W/V =",
                       format(c(W(mod1) / V(mod1),
                                W(mod2) / V(mod2)))))
legend("bottomleft", legend = leg,
       col = c("darkgrey", "black", "black"),
       lty = c("solid", "longdash", "dotdash"),
       pch = c(1, NA, NA), bty = "n")

###
mod0 <- dlmModPoly(order = 1, dV = 15100, dW = 1468)
X <- ts(matrix(mod0$W, nc = 1, nr = length(Nile)),
        start = start(Nile))
window(X, 1898, 1899) <- 12 * mod0$W
modDam <- mod0
modDam$X <- X
modDam$JW <- matrix(1, 1, 1)
damFilt <- dlmFilter(Nile, modDam)
mod0Filt <- dlmFilter(Nile, mod0)
a <- window(cbind(Nile, mod0Filt$f, damFilt$f),
            start = 1880, end = 1920)
plot(a[, 1], type = 'o', col = "darkgrey",
     xlab = "", ylab = "Level")
lines(a[, 2], lty = "longdash")
lines(a[, 3], lty = "dotdash")
abline(v = 1898, lty = 2)
leg <- c("data", paste("one-step-ahead forecast -",
                       c("mod0", "modDam")))
legend("bottomleft", legend = leg,
       col = c("darkgrey", "black", "black"),
       lty = c("solid", "longdash", "dotdash"),
       pch = c(1, NA, NA), bty = "n")

###
set.seed(1)
expdFore <- dlmForecast(expdFilt, nAhead = 12, sampleNew = 10)
plot(window(expd, start = c(1964, 1)), type = 'o',
     xlim = c(1964, 1971), ylim = c(350, 850),
     xlab = "", ylab = "Expenditures")
names(expdFore)
attach(expdFore)
invisible(lapply(newObs, function(x)
                 lines(x, col = "darkgrey",
                       type = 'o', pch = 4)))
lines(f, type = 'o', lwd = 2, pch = 16)
abline(v = mean(c(time(f)[1], time(expd)[length(expd)])),
       lty = "dashed")
detach()

###
qqnorm(residuals(damFilt, sd = FALSE))
qqline(residuals(damFilt, sd = FALSE))
tsdiag(damFilt)
