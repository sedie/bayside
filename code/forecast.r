args <- commandArgs(trailingOnly=TRUE)
ftime <- as.numeric(args[1])

library(rstan)
library(plyr)
library(parallel)
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
    outs <- extract(zips, permuted = TRUE) # sample the posterior
    p <- data$P # num of groups
    initial <- sapply(1:nrow(data$counts), function(i) data$counts[ i, data$end[i]] )

    # regression
    coef0 <- apply(outs$coef[, , 1], 2, function(x) sample(x, 1))
    coef1 <- apply(outs$coef[, , 2], 2, function(x) sample(x, 1))

    # acp
    alp <- apply(outs$alpha, 2, function(x) sample(x, 1))  # for each prov
    bet <- apply(outs$beta, 2, function(x) sample(x, 1))  # for each prov

    # markov
    gam <- apply(outs$gamma, 2, function(x) sample(x, 1))  # for each prov
    eta <- apply(outs$eta, 2, function(x) sample(x, 1))  # for each prov

    # by province
    out <- list()
    for(ii in seq(p)) {
        all_toff <- rev(rev(data$off[ii, ])[seq(data$end[ii])])
        # generate offset segment by sampling past decade of offsets
        toff <- sample(all_toff[ (length(all_toff) - ftime ) : length(all_toff) ], ftime, replace=TRUE)

        # by time point
        mu <- c()
        theta <- c()
        oo <- data$counts[ii, data$end[ii]]
        co <- c() # empty expected count vec

        for(jj in 1:ftime) {  # I've changed the time seq to extend the series by ftime years
          # jj = 2
          if(jj == 1) { # the starting point is the number of species described at end of time series.
            mu[jj] <- oo
            theta[jj] <- 0
            co[jj] <- oo[jj]
          } else { # if not the starting point
            mu[jj] <- exp(coef0[ii] + coef1[ii] * jj) +
            alp[ii] * oo[jj - 1] + bet[ii] * mu[jj - 1]

            val <- (oo[jj - 1] == 0) * 1
            theta[jj] <- (val * gam[ii]) + ((1 - val) * eta[ii])

            # need to add in offset (# number of papers per year)
            co[jj] <- rZIP(1, mu = (toff[jj] + 1) * mu[jj], 
                           sigma = theta[jj])
            # update observed count to the new expected count - seems necessary if forecasting
            oo[jj] <- co[jj]
          }
        }
        out[[ii]] <- co
    }
    out
}

# simulate the forecast
forecast <- mclapply(1:1000, mc.cores=4, function(ii) {
   post.forecast(data = data, ftime=ftime, model = zips) 
})
save(forecast, file="data/dump/forecast.data")