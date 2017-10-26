library(ggplot2)
library(forecast)
library(dlm)
library(zoo)
library("quantmod")
library(tseries)

amd <- getSymbols("AMD", from="2017-01-01", to="2017-10-01", auto.assign = FALSE)$AMD.Adjusted
autoplot(amd)
nvda <- getSymbols("NVDA", from="2017-01-01", to="2017-10-01", auto.assign = FALSE)$NVDA.Adjusted
autoplot(nvda)
smh <- getSymbols("SMH", from="2017-01-01", to="2017-10-01", auto.assign = FALSE)$SMH.Adjusted
autoplot(smh)
mu <- getSymbols("MU", from="2017-01-01", to="2017-10-01", auto.assign = FALSE)$MU.Adjusted
autoplot(mu)
ibm <- getSymbols("IBM", from="2017-01-01", to="2017-10-01", auto.assign = FALSE)$IBM.Adjusted
autoplot(ibm)

df <- data.frame(index(amd), amd, nvda, mu, smh, ibm)
(tss <- read.zoo(df))
autoplot(tss) + facet_grid(Series ~ ., scales = "free_y")

lm(nvda ~ mu) %>% summary()
lm(nvda ~ smh) %>% summary()
lm(nvda ~ amd) %>% summary()
lm(nvda ~ ibm) %>% summary()
lm(nvda ~ ibm) %>% residuals() %>% qqnorm()

fit <- lm(nvda ~ ibm)
adf.test(fit$residuals)

### dlm regression ###

predictor <- ibm
build <- function(u) {
  dlmModReg(predictor, dV = exp(u[1]), dW = exp(u[2 : 3]))
}
outMLE <- dlmMLE(nvda, parm = rep(0, 3), build)
outMLE$convergence
exp(outMLE$par)
outMLE$value
mod <- build(outMLE$par)
mod$FF
mod$GG
mod$W
mod$V

outF <-dlmFilter(nvda, mod)
outF$m[1 + length(nvda), ]
outS <- dlmSmooth(nvda, mod)


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



