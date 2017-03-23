args <- commandArgs(trailingOnly=TRUE)
ftime <- as.numeric(args[1])

library(rstan)
library(plyr)
library(gamlss.dist)

set.seed(420)

# initial data
files <- dir('data/dump', pattern = 'count_info.data.R',
             full.names = TRUE)
data <- read_rdump(files)

# zero inflated fits
load("data/dump/fit.data") # loads as fit
zips <- fit # reassign to zips
rm(fit) # remove from memory

post.forecast <- function(data, ftime, model) {
    outs <- extract(model, permuted = TRUE) # posterior
    p <- data$P # number of groups

    # regression
    coef0 <- apply(outs$coef[, , 1], 2, function(x) sample(x, 1))
    coef1 <- apply(outs$coef[, , 2], 2, function(x) sample(x, 1))

    # acp
    alp <- apply(outs$alpha, 2, function(x) sample(x, 1))  # for each prov
    bet <- apply(outs$beta, 2, function(x) sample(x, 1))  # for each prov

    # markov
    gam <- apply(outs$gamma, 2, function(x) sample(x, 1))  # for each prov
    eta <- apply(outs$eta, 2, function(x) sample(x, 1))  # for each prov

    # t=1
    phi <- apply(outs$phi, 2, function(x) sample(x, 1))  # for each prov

    # initial count
    initial <- sapply(seq(p), function(pp) data$counts[ pp, data$str[pp]])

    # by province
    out <- list()
    for(ii in seq(p)) {
        all_toff <- data$off[ii, ][data$str[ii]:data$end[ii]] # offset starts at first naming year
        # generate offset segment by sampling past decade of offsets
        toff <- sample(all_toff[ (length(all_toff) - ftime ) : length(all_toff) ], 
            ftime, replace=TRUE)

        mu <- c()
        theta <- c()
        oo <- data$counts[ii, data$end[ii]]
        co <- c() # empty expected count vec
        # by time point
        for(jj in 1:ftime) {  # time seq extends the series by ftime years
          if(jj == 1) { # the starting point is the number of species described at end of time series.
            mu[jj] <- oo
            theta[jj] <- 0
            co[jj] <- oo[jj]
          } else { # if not the starting point
            mu[jj] <- exp(coef0[ii] + coef1[ii] * jj) +
            alp[ii] * oo[jj - 1] + bet[ii] * mu[jj - 1]

            val <- (oo[jj - 1] == 0) * 1
            theta[jj] <- (val * gam[ii]) + ((1 - val) * eta[ii])

            co[jj] <- rZIP(1, mu = (toff[jj] + 1) * mu[jj], 
                           sigma = theta[jj])
            oo[jj] <- co[jj] # current count becomes next expected count
          }
        }
        out[[ii]] <- co
    }
    out
}

# simulate the forecast
forecast <- lapply(1:1000, function(ii) {
   post.forecast(data=data, ftime=ftime, model=zips) 
})
save(forecast, file="data/dump/forecast.data")