### R code for Chapter 5
require(dlm)


###
### Basic PF: A simulated local level model
###

### Generate data
mod <- dlmModPoly(1, dV = 2, dW = 1, m0 = 10, C0 = 9)
n <- 100
set.seed(23)
simData <- dlmForecast(mod = mod, nAhead = n, sampleNew = 1)
y <- simData$newObs[[1]]

### Basic Particle Filter - optimal importance density
N <- 1000
N_0 <- N / 2
pfOut <- matrix(NA_real_, n + 1, N)
wt <- matrix(NA_real_, n + 1, N)
importanceSd <- sqrt(drop(W(mod) - W(mod)^2 / (W(mod) + V(mod))))
predSd <- sqrt(drop(W(mod) + V(mod)))
## Initialize sampling from the prior
pfOut[1, ] <- rnorm(N, mean = m0(mod), sd = sqrt(C0(mod)))
wt[1, ] <- rep(1/N, N)
for (it in 2 : (n + 1))
{
    ## generate particles
    means <- pfOut[it - 1, ] + W(mod) *
        (y[it - 1] - pfOut[it - 1, ]) / (W(mod) + V(mod))
    pfOut[it, ] <- rnorm(N, mean = means, sd = importanceSd)
    ## update the weights
    wt[it, ] <- dnorm(y[it - 1], mean = pfOut[it - 1, ], sd = predSd) *
        wt[it - 1, ]
    wt[it, ] <- wt[it, ] / sum(wt[it, ])
    ## resample, if needed
    N.eff <- 1 / crossprod(wt[it, ])
    if ( N.eff < N_0 )
    {
        ## multinomial resampling
        index <- sample(N, N, replace = TRUE, prob = wt[it, ])
        pfOut[it, ] <- pfOut[it, index]
        wt[it, ] <- 1 / N
    }
}

### Plot results
### Compare exact filtering distribution with PF approximation
modFilt <- dlmFilter(y, mod)
thetaHatKF <- modFilt$m[-1]
sdKF <- with(modFilt, sqrt(unlist(dlmSvd2var(U.C, D.C))))[-1]
pfOut <- pfOut[-1, ]
wt <- wt[-1, ]
thetaHatPF <- sapply(1 : n, function(i) weighted.mean(pfOut[i, ], wt[i, ])) 
sdPF <- sapply(1 : n, function(i)
               sqrt(weighted.mean((pfOut[i, ] - thetaHatPF[i])^2, wt[i, ])))
par(mfrow = c(2, 1))
plot.ts(cbind(thetaHatKF, thetaHatPF),
        plot.type = "s", lty = c("dotted", "longdash"),
        xlab = "", ylab = expression(m[t]))
legend("topleft", c("Kalman", "Particle"),
       lty = c("dotted", "longdash"), bty = "n")
plot.ts(cbind(sdKF, sdPF), plot.type = "s",
        lty = c("dotted", "longdash"), xlab = "",
        ylab = expression(sqrt(C[t])))
legend("topright", c("Kalman", "Particle"),
       lty = c("dotted", "longdash"), bty = "n")

###
### PF with unknown parameters: Liu and West
### (use same data as before, but assume unknown V and W)
###

N <- 10000
a <- 0.975
set.seed(4521)
pfOutTheta <- matrix(NA_real_, n + 1, N)
pfOutV <- matrix(NA_real_, n + 1, N)
pfOutW <- matrix(NA_real_, n + 1, N)
wt <- matrix(NA_real_, n + 1, N)
## Initialize, sampling from the prior
pfOutTheta[1, ] <- rnorm(N, mean = m0(mod),
                         sd = sqrt(C0(mod)))
pfOutV[1, ] <-  runif(N, 0, 10)
pfOutW[1, ] <-  runif(N, 0, 10)
wt[1, ] <- rep(1/N, N)
for (it in 2 : (n + 1))
{
    ## compute means and variances of the particle
    ## cloud for V and W
    meanV <- weighted.mean(pfOutV[it - 1, ], wt[it - 1,])
    meanW <- weighted.mean(pfOutW[it - 1, ], wt[it - 1,])
    varV <- weighted.mean((pfOutV[it - 1, ] - meanV)^2,
                          wt[it - 1,])
    varW <- weighted.mean((pfOutW[it - 1, ] - meanW)^2,
                          wt[it - 1,])
    ## compute the parameters of Gamma kernels
    muV <- a * pfOutV[it - 1, ] + (1 - a) * meanV
    sigma2V <- (1 - a^2) * varV
    alphaV <- muV^2 / sigma2V
    betaV <- muV / sigma2V
    muW <- a * pfOutW[it - 1, ] + (1 - a) * meanW
    sigma2W <- (1 - a^2) * varW
    alphaW <- muW^2 / sigma2W
    betaW <- muW / sigma2W
    ## draw the auxiliary indicator variables
    probs <- wt[it - 1,] * dnorm(y[it - 1], sd = sqrt(muV),
                                 mean = pfOutTheta[it - 1, ])
    auxInd <- sample(N, N, replace = TRUE, prob = probs)
    ## draw the variances V and W
    pfOutV[it, ] <- rgamma(N, shape = alphaV[auxInd],
                           rate = betaV[auxInd])
    pfOutW[it, ] <- rgamma(N, shape = alphaW[auxInd],
                           rate = betaW[auxInd])
    ## draw the state theta
    pfOutTheta[it, ] <- rnorm(N, mean =
                              pfOutTheta[it - 1, auxInd],
                              sd = sqrt(pfOutW[it, ]))
    ## compute the weights
    wt[it, ] <- exp(dnorm(y[it - 1], mean = pfOutTheta[it, ],
                          sd = sqrt(pfOutV[it, ]),
                          log = TRUE) -
                    dnorm(y[it - 1], mean =
                          pfOutTheta[it - 1, auxInd],
                          sd = sqrt(muV[auxInd]),
                          log = TRUE))
    wt[it, ] <- wt[it, ] / sum(wt[it, ])
}

## Look at 'filtering' distributions: mean, and 50% prob interval
pfOutV <- pfOutV[-1, ]
pfOutW <- pfOutW[-1, ]
wt <- wt[-1, ]

meanV <- mapply(weighted.mean, split(pfOutV, 1 : n),
                split(wt, 1 : n))
meanW <- mapply(weighted.mean, split(pfOutW, 1 : n),
                split(wt, 1 : n))

### require function "weighted.quantile"
quantV <- t(mapply(weighted.quantile,
                   split(pfOutV, 1 : n), split(wt, 1 : n),
                   MoreArgs = list(probs = c(0.25, 0.75))))
quantW <- t(mapply(weighted.quantile,
                   split(pfOutW, 1 : n), split(wt, 1 : n),
                   MoreArgs = list(probs = c(0.25, 0.75))))
par(mfrow = c(2, 1))
plot.ts(cbind(meanV, quantV), plot.type = "s",
        lty = c("solid", rep("longdash", 2)),
        xlab = "", ylab = "V")
abline(h = V(mod), lty = "dotted")
legend("topright", c("posterior mean",
                     "50% probability interval", "true"),
       lty = c("solid", "longdash", "dotted"), bty = "n")
plot.ts(cbind(meanW, quantW), plot.type = "s",
        lty = c("solid", rep("longdash", 2)),
        xlab = "", ylab = "W")
abline(h = W(mod), lty = "dotted")
legend("topright", c("posterior mean",
                     "50% probability interval", "true"),
       lty = c("solid", "longdash", "dotted"), bty = "n")

