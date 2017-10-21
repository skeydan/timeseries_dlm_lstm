library(ggplot2)
library(forecast)
library(dlm)
library(zoo)
library("quantmod")

nvda <- getSymbols("NVDA", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$NVDA.Adjusted
autoplot(nvda)

# see https://stats.stackexchange.com/questions/214581/why-is-my-kalman-filter-trusting-so-much-my-observations!!

(index <- index(nvda))
nvda <- as.ts(nvda)

StructTS(nvda)

reldiff <- function(vec) diff(vec)/vec[-length(vec)]
nvda <- reldiff(nvda)

(log_var_diff <- log(var(diff(nvda))))


# local level

build <- function(params) {
  dlmModPoly(order = 1, dV = exp(params[1]), dW = exp(params[2]))
}
initial_params <- rep(log_var_diff, 2)
ll_fit <- dlmMLE(nvda, initial_params, build)
ll_fit$convergence

model <- dlmModPoly(order = 1, dV = exp(ll_fit$par[1]), dW = exp(ll_fit$par[2]))
#model <- dlmModPoly(order = 1, dV = 100, dW = 0.5)

unlist(model)

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(nvda, col = "darkgrey", xlim = c(0, 500))
filt <- dlmFilter(nvda, model)
lines(dropFirst(filt$m), lty = "longdash", col="red")
smooth <- dlmSmooth(nvda, model)
lines(smooth$s, lty = "dotdash", col="green")
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue", lwd=.4)

### Diagnostic plots
res <- residuals(filt, sd=FALSE)
plot(res,type='h'); abline(h=0)
acf(res)
qqnorm(res); qqline(res)


### filtering & smoothing estimates - intervals
par(mfrow = c(2, 2), mar = c(3, 4, 1, 1) + 0.1, cex = 0.5)
plot(nvda, xlab = "", ylab = "(a)", col = "darkgrey")
lines(dropFirst(filt$m), lty = "longdash")
sdFilt <- with(filt, ts(sqrt(unlist(dlmSvd2var(U.C, D.C))), start = start(y)))
lns <- cbind(dropFirst(filt$m), qnorm(0.95, sd = sdFilt[-1]) %o% c(-1, 1) + filt$m[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash")
plot(dropFirst(sdFilt), xlab = "", ylab = "(b)")
plot(nvda, xlab = "", ylab = "(c)", col = "darkgrey")
lines(dropFirst(smooth$s), lty = "longdash")
sdSmooth <- with(smooth, ts(sqrt(unlist(dlmSvd2var(U.S, D.S))),start = start(nvda)))
lns <- cbind(dropFirst(smooth$s), qnorm(0.95, sd = sdSmooth[-1]) %o% c(-1, 1) + smooth$s[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash") 
plot(dropFirst(sdSmooth), xlab = "", ylab = "(d)")

### Formal diagnostic tests
shapiro.test(res)
Box.test(res, lag=20, type="Ljung")
sapply(1 : 20, function(i)
  Box.test(res, lag = i, type = "Ljung-Box")$p.value)


# linear trend

build <- function(params) {
  dlmModPoly(order = 2, dV = exp(params[1]), dW = exp(params[2:3]))
}

initial_params <- rep(log_var_diff, 3)
lg_fit <- dlmMLE(nvda, initial_params, build)
lg_fit$convergence

model <- dlmModPoly(order = 2, dV = exp(lg_fit$par[1]), dW = exp(c(lg_fit$par[2], lg_fit$par[3])))
model$FF
model$GG
model$W
model$V

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(nvda, type='o', col = "darkgrey", xlim = c(0, 500))
filt <- dlmFilter(nvda, model)
lines(dropFirst(filt$m[ ,1]), lty = "longdash", col="red") # column 1 of state
smooth <- dlmSmooth(nvda, model)
lines(smooth$s[ ,1], lty = "dotdash", col="green") # column 1 of state
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue")



### model comparisons
mod1 <- dlmModPoly(order = 1, dV = 0.4, dW = 6.8)
mod2 <- dlmModPoly(order = 2, dV = 0.5, dW = c(6.5, 0.00000001) )
filt1 <- dlmFilter(nvda, mod1)
fut1 <- dlmForecast(filt, n=12)
filt2 <- dlmFilter(nvda, mod2)
fut2 <- dlmForecast(filt2, n=12)

mean(abs(filt1$f - nvda)) # MAD
mean(abs(filt2$f - nvda))
mean((filt1$f - nvda)^2) # MSE
mean((filt2$f - nvda)^2)
mean(abs(filt1$f - nvda) / nvda) # MAPE
mean(abs(filt2$f - nvda) / nvda)





###

# ewgs <- getSymbols("EWGS", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$EWGS.Adjusted
# autoplot(ewgs)
# 
# eurl <- getSymbols("EURL", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$EURL.Adjusted
# autoplot(eurl)
# 
# ewo <- getSymbols("EWO", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$EWO.Adjusted
# autoplot(ewo)
# 
# smez <- getSymbols("SMEZ", from="2016-01-01", to="2017-10-01", auto.assign = FALSE)$SMEZ.Adjusted
# autoplot(smez)


