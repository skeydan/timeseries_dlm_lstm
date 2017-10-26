source("get_olympics_data.R")
library(ggplot2)
library(forecast)
library(dlm)

ggplot(male400, aes(x = year, y = seconds)) + geom_line() +
  geom_point(data = male400, size=1) +
  ggtitle("Men's 400m Olympic winning times 1986-2016") +
  xlab("year") + ylab("seconds") 

StructTS(male400$seconds)

### linear regression approach ###
fit <- lm(seconds ~ year, male400_1996)
summary(fit)
acf(fit$residuals)
qqnorm(fit$residuals) #!!!

df <- male400_1996 %>% bind_cols(is_pred = factor(rep(0, nrow(male400_1996)), levels = c(0,1)))
df <-  df %>% bind_rows(data.frame(year = seq(2000, 2016, by=4), is_pred = factor(rep(1,5), levels = c(0,1))))
preds <- fit %>% predict(newdata = df, interval = "prediction")
df[25:29,2] <- preds[25:29,1]
ggplot(df, aes(x = year, y=seconds)) + geom_point(aes(color = is_pred, shape = is_pred), size = 1) + 
  geom_ribbon(aes(ymin = preds[ , 2], ymax = preds[ , 3]), alpha = 0.2) +
  ggtitle("Men's 400m Olympic winning times 1986-1996 and predictions for 2000-2016", subtitle = "Prediction intervals from least squares") +
  xlab("year") + ylab("seconds") 


### NO arima because this is not a regular timeseries!!!


### dlm regression ###

build <- function(u) {
  dlmModReg(male400_1996$year, dV = exp(u[1]), dW = exp(u[2 : 3]))
}
outMLE <- dlmMLE(male400_1996$seconds, parm = rep(0, 3), build)
outMLE$convergence
exp(outMLE$par)
outMLE$value
mod <- build(outMLE$par)
mod$FF
mod$GG
mod$W

outF <-dlmFilter(male400_1996$seconds, mod)
outF$m[1 + length(male400_1996$seconds), ]
outS <- dlmSmooth(male100_1996$seconds, mod)


require(zoo)
# development of coefficients
plot(as.zoo(dropFirst(outF$m)), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5)
plot(as.zoo(dropFirst(outS$s)), main = "", mar = c(0, 2.1, 0, 1.1),
     oma = c(2.1,0,.1,.1), cex.axis = 0.5)

# forecasts?
dlmForecast(outF, n=5)
# dlmForecast only works with constant models
# Constant refers to the transition, observation matrices, not the state itself.
# In dlm terms, you have a constant model if you are not using the "J" components like JGG or JFF (if I remember correctly).

mod$JFF
mod$JGG
mod$JV
mod$JW
mod$X
