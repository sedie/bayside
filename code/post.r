# run posterior simulation from model fit
library(rstan)
library(gamlss.dist)
library(parallel)

set.seed(420)

# initial data
files <- dir('data/dump', pattern = 'count_info.data.R',
             full.names = TRUE)
data <- read_rdump(files)

# zero inflated fits
load("data/dump/fit.data") # loads as fit
zips <- fit # reassign to zips
rm(fit) # remove from memory

# posterior predictive simulations for checking model fit
posterior.sim <- function(data, model, over = FALSE) {
  sims <- list()
  outs <- extract(model, permuted = TRUE)
  p <- data$P
  initial <- data$counts[, 1]

  # regression
  coef0 <- apply(outs$coef[, , 1], 2, function(x) sample(x, 1))
  coef1 <- apply(outs$coef[, , 2], 2, function(x) sample(x, 1))

  # acp
  alp <- apply(outs$alpha, 2, function(x) sample(x, 1))  # for each prov
  bet <- apply(outs$beta, 2, function(x) sample(x, 1))  # for each prov

  # markov
  gam <- apply(outs$gamma, 2, function(x) sample(x, 1))  # for each prov
  eta <- apply(outs$eta, 2, function(x) sample(x, 1))  # for each prov

  if(over) phi <- apply(outs$phi, 2, function(x) sample(x, 1))

  # by province
  for(ii in seq(p)) {
    # get the offset segment
    toff <- rev(rev(data$off[ii, ])[seq(data$end[ii])])

    # by time point
    mu <- c()
    theta <- c()
    oo <- c()  # counts over time
    for(jj in seq(data$end[ii])) {
      if(jj == 1) {
        mu[jj] <- exp(coef0[ii] + coef1[ii] * jj)
        theta[jj] <- 0

        oo[jj] <- initial[ii]
      } else {
        mu[jj] <- exp(coef0[ii] + coef1[ii] * jj) +
        alp[ii] * oo[jj - 1] + bet[ii] * mu[jj - 1]

        val <- (oo[jj - 1] == 0) * 1
        theta[jj] <- (val * gam[ii]) + ((1 - val) * eta[ii])

        if(!over) {
          # need to add in offset (# number of papers per year)
          oo[jj] <- rZIP(1, mu = (toff[jj] + 1) * mu[jj], 
                         sigma = theta[jj])
        } else if(over) {
          if(runif(1, min = 0, max = 1) > theta[jj]) {
            oo[jj] <- rnbinom(1, mu = (toff[jj] + 1) * mu[jj], 
                              size = phi[ii])
          } else {
            oo[jj] <- 0
          }

        }
      }
    }
    sims[[ii]] <- oo
  }
  sims
}

# post sim
allsim <- mclapply(1:1000, mc.cores=4, function(ii) {
    posterior.sim(data = data, model = zips, over = FALSE)
} )
save(allsim, file="data/dump/post.data")