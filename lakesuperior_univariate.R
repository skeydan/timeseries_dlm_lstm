library(ggplot2)
library(forecast)
library(dlm)

lakeSup <- ts(read.table("Datasets/lakeSuperior.dat", skip = 3,
                         colClasses = "numeric")[, 2], start = 1900)

autoplot(lakeSup)
StructTS(lakeSup)

(log_var_diff <- log(var(diff(lakeSup))))


# local level

build <- function(params) {
  dlmModPoly(order = 1, dV = exp(params[1]), dW = exp(params[2]))
}
initial_params <- rep(log_var_diff, 2)
ll_fit <- dlmMLE(lakeSup, initial_params, build)
ll_fit$convergence

model <- dlmModPoly(order = 1, dV = exp(ll_fit$par[1]), dW = exp(ll_fit$par[2]))
unlist(model)

### filter - smooth - forecast
par(mfrow = c(1,1))
plot(lakeSup, col = "darkgrey", xlim = c(1900, 2000))
filt <- dlmFilter(lakeSup, model)
lines(dropFirst(filt$m), lty = "longdash", col="red")
smooth <- dlmSmooth(lakeSup, model)
lines(smooth$s, lty = "dotdash", col="green")
fc <- dlmForecast(filt, nAhead = 12)
lines(fc$f, type = "o", col="blue", lwd=.4)

### Diagnostic plots
res <- residuals(filt, sd=FALSE)
plot(res,type='h'); abline(h=0)
acf(res)
qqnorm(res); qqline(res)


