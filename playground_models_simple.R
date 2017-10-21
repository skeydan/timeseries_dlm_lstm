library(dlm)
library(ggplot2)
library(forecast)

#y <- as.ts(c(1:10))
y <- AirPassengers
#y <- window(AirPassengers, start = 1958)

(log_var_y_diff <- log(var(diff(y))))

StructTS(y)

############################## local level model #####################################


build <- function(params) {
  dlmModPoly(order = 1, dV = exp(params[1]), dW = exp(params[2]))
}
initial_params <- rep(log_var_y_diff, 2)
ll_fit <- dlmMLE(y, initial_params, build)
ll_fit$convergence

model <- dlmModPoly(order = 1, dV = exp(ll_fit$par[1]), dW = exp(ll_fit$par[2]))
unlist(model)

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(y, type='o', col = "darkgrey", xlim = c(1958, 1962))
filt <- dlmFilter(y, model)
lines(dropFirst(filt$m), lty = "longdash", col="red")
smooth <- dlmSmooth(y, model)
lines(smooth$s, lty = "dotdash", col="green")
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue")

### Diagnostic plots
res <- residuals(filt, sd=FALSE)
plot(res,type='h'); abline(h=0)
acf(res)
qqnorm(res); qqline(res)


### filtering & smoothing estimates - intervals
par(mfrow = c(2, 2), mar = c(3, 4, 1, 1) + 0.1, cex = 0.5)
plot(y, xlab = "", ylab = "(a)", col = "darkgrey")
lines(dropFirst(filt$m), lty = "longdash")
sdFilt <- with(filt, ts(sqrt(unlist(dlmSvd2var(U.C, D.C))), start = start(y)))
lns <- cbind(dropFirst(filt$m), qnorm(0.95, sd = sdFilt[-1]) %o% c(-1, 1) + filt$m[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash")
plot(dropFirst(sdFilt), xlab = "", ylab = "(b)")
plot(y, xlab = "", ylab = "(c)", col = "darkgrey")
lines(dropFirst(smooth$s), lty = "longdash")
sdSmooth <- with(smooth, ts(sqrt(unlist(dlmSvd2var(U.S, D.S))),start = start(y)))
lns <- cbind(dropFirst(smooth$s), qnorm(0.95, sd = sdSmooth[-1]) %o% c(-1, 1) + smooth$s[-1]) 
for (i in 1 : 2) lines(lns[, i + 1], lty = "dotdash") 
plot(dropFirst(sdSmooth), xlab = "", ylab = "(d)")

### Formal diagnostic tests
shapiro.test(res)
Box.test(res, lag=20, type="Ljung")
sapply(1 : 20, function(i)
  Box.test(res, lag = i, type = "Ljung-Box")$p.value)


####################### linear growth model ##################################

build <- function(params) {
  dlmModPoly(order = 2, dV = exp(params[1]), dW = exp(params[2:3]))
}
initial_params <- rep(log_var_y_diff, 3)
lg_fit <- dlmMLE(y, initial_params, build)
lg_fit$convergence

model <- dlmModPoly(order = 2, dV = exp(lg_fit$par[1]), dW = exp(c(lg_fit$par[2], lg_fit$par[3])))
model$FF
model$GG
model$W
model$V

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(y, type='o', col = "darkgrey", xlim = c(1958, 1962))
filt <- dlmFilter(y, model)
lines(dropFirst(filt$m[ ,1]), lty = "longdash", col="red") # column 1 of state
smooth <- dlmSmooth(y, model)
lines(smooth$s[ ,1], lty = "dotdash", col="green") # column 1 of state
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue")


### model comparisons
mod1 <- dlmModPoly(order = 1, dV = ll_fit$par[1], dW = c(ll_fit$par[2]))
mod2 <- dlmModPoly(order = 2, dV = lg_fit$par[1], dW = c(lg_fit$par[2], lg_fit$par[3]))
filt1 <- dlmFilter(y, mod1)
fut1 <- dlmForecast(filt, n=12)
filt2 <- dlmFilter(y, mod2)
fut2 <- dlmForecast(filt2, n=12)

mean(abs(filt1$f - y)) # MAD
mean(abs(filt2$f - y))
mean((filt1$f - y)^2) # MSE
mean((filt2$f - y)^2)
mean(abs(filt1$f - y) / y) # MAPE
mean(abs(filt2$f - y) / y)

# Theil's u
# if U < 1 model produces better forecasts than trivial "no change" model
sqrt(sum((filt1$f - y)[-(1:5)]^2) / sum(diff(y[-(1:4)])^2))
sqrt(sum((filt2$f - y)[-(1:5)]^2) / sum(diff(y[-(1:4)])^2))

####################### Seasonal models ##################################

#build <- function(params) {
#  dlmModSeas(frequency = 12, dV = params[1], dW = params[2:12])
#}
build <- function(params) {
  dlmModSeas(frequency = 4, dV = exp(params[1]), dW = exp(params[2:4]))
}

#initial_params <- rep(var_y_diff, 12)
initial_params <- rep(log_var_y_diff, 4)

ss_fit <- dlmMLE(y, initial_params, build)
ss_fit$convergence

model <-  dlmModSeas(frequency = 12, dV = exp(lg_fit$par[1]), dW = exp(lg_fit$par[2:12]))
model$FF
model$GG
model$W
model$V

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(y, type='o', col = "darkgrey", xlim = c(1958, 1962))
filt <- dlmFilter(y, model, debug = TRUE)


####################### ARIMA models ##################################

auto.arima(y)
build <- function(params) {
  dlmModARMA(ar = params[1:2], ma = params[3], sigma2 = params[4], dV = 0.001)
}

initial_params <- c(0.6, 0.2, -1, 132)

arma_fit <- dlmMLE(y, initial_params, build)
arma_fit$convergence
unlist(build(arma_fit$par))

model <-  dlmModARMA(ar = arma_fit$par[1:2], ma = arma_fit$par[3], sigma2 = exp(arma_fit$par[4]), dV = 0.001)
model$FF
model$GG
model$W
model$V

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(y, type='o', col = "darkgrey", xlim = c(1958, 1962))
filt <- dlmFilter(y, model, debug = TRUE)
lines(dropFirst(filt$m[ ,1]), lty = "longdash", col="red") # column 1 of state
smooth <- dlmSmooth(y, model)
lines(smooth$s[ ,1], lty = "dotdash", col="green") # column 1 of state
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue")


