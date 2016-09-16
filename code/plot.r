# visualize results from model fit and forecasting
args <- commandArgs(trailingOnly=TRUE)
OFFSET <- as.logical(args[1])
# OFFSET <- TRUE

library(dplyr)
library(plyr)
library(rstan)
library(ggplot2)

#
#    LOAD DATA
#

# original data
dataRAW <- read.csv("data/data.csv")

# initial data
files <- dir('data/dump', pattern = 'count_info.data.R',
             full.names = TRUE)
data <- read_rdump(files)

# zero inflated fits
load("data/dump/fit.data") # loads as fit
zips <- fit # reassign to zips
rm(fit) # remove from memory

# posterior samples
load("data/dump/post.data") # loads as allsim

# forecast data
load("data/dump/forecast.data") # loads as forecast

# map model indeces to original variables
mapping <- data.frame( 
    groupname=as.character(unique(dataRAW$group)),
    group=as.numeric(unique(dataRAW$group)))

#
#   END \\ LOAD
#


#
#   PREPARE CUMULATIVE SERIES
#

# cumulative series for observed data
cumm <- lapply( seq(nrow(data$counts)), function(ii) { # each prov
    data.frame(
            index=data$str[ii]:data$end[ii],
             value=data$counts[ii, 1:data$end[ii]],
             cml_value=cumsum(data$counts[ii, 1:data$end[ii]]),
             off=data$off[ii, ][ data$str[ii]:data$end[ii] ],
             group=ii
    )
}) %>% rbind.fill
cumm$sim <- 0 # set sim to 0 to mark observed data

# cumulative series for sampled data
cummsim <- lapply( seq(length(allsim)), function(jj) { # each sim
    lapply( seq(length(data$end)), function(ii) { # each prov
        data.frame(
            index=ncol(data$counts) - rev(seq(data$end[ii])), 
            value=allsim[[jj]][[ii]],
            cml_value=cumsum(allsim[[jj]][[ii]]), group=ii, sim=jj,
            off=data$off[ii, ][ data$str[ii]:data$end[ii] ],
        )
    } ) %>% rbind.fill
}) %>% rbind.fill

both <- rbind(cumm, cummsim) # combine both the off and nooff data
Z <- split(both, both$group) # split by group

#
#   END \\ PREPARE CUMULATIVE SERIES
#



#
#   SUMMARIZE MODEL RESULTS
#


# for each group, determine the mean number of species per simulation
ex.med <- lapply( 1:length(allsim[[1]]), function(ii) {
    median( unlist( lapply( sapply(allsim, "[", ii), sum) ) )
} ) %>% do.call(rbind, .)
# and the CI
ex.CI <- lapply( 1:length(allsim[[1]]), function(ii){
    quantile( unlist( lapply( sapply(allsim, "[", ii), sum) ), probs=c(0.1,0.9) )
} ) %>% do.call(rbind, . )
# and the observed count
obs.count <- lapply(Z, function(z) {
    max( z[ z$sim==0, "cml_value"] )
}) %>% do.call(rbind, .)

# outs$coef[2] is the coefficient that estimates the long-term trend in \\
# description rate, it's in log units so need to exponentiate it \\
# where 0 equals means no trend, positive means increasing description \\
#  through time and negative is decreasing...
# `coef` is an array where \\
# [posterior sample [1:4000], province number[1:16], coefficient[1:2] ]
coef <- rstan::extract(zips, par="coef")[[1]]
# pull out the second coefficient for each province across each 
# posterior sample
cf2 <- apply(coef, 1, function(i) i[ , 2] )
median_cf2 <- data.frame(provid=1:nrow(cf2),
    slowdown=apply(cf2, 1, median) )
CI80_cf2 <- data.frame(provid=1:nrow(cf2), 
    slowdown=t( apply(cf2, 1, quantile, probs=c(0.1, 0.9))) )

# compile results into table
results <- data.frame(
    group=mapping[ match(as.numeric(rownames(obs.count)), mapping[, 2]), 2 ] , 
    groupname=mapping[ match(as.numeric(rownames(obs.count)), mapping[, 2]), 1 ] , 
    observed_species=obs.count , 
    expected_median=ex.med, 
    expected_CI_lower=ex.CI[,1], 
    expected_CI_higher=ex.CI[,2],
    slowdown=round(signif(median_cf2$slowdown, 2 ), 4 ),
    slowdown_CI_lower=round(signif(CI80_cf2[,2], 2), 4),
    slowdown_CI_higher=round(signif(CI80_cf2[,3], 2), 4)
)


#   Long-term regression lines

calc.lambda <- function(x, b0, b1) {
  exp(b0 + b1 * x)
}

# pull out ends of time sequences
end <- data$end
# pull out the second coefficient for each province across each posterior sample
prov.cf1 <- apply(coef, 1, function(i) i[ ,1])
prov.cf2 <- apply(coef, 1, function(i) i[ ,2])
# for each group
lambda <- lapply(1:length(end), function(ii) {   
    # ii <- 1
    cf1 <- prov.cf1[ii, ]
    cf2 <- prov.cf2[ii, ]
    # for each coefficient pair
    cPair <- lapply(1:length(cf1), function(kk){
        # kk <- 1
        b0 <- cf1[kk]
        b1 <- cf2[kk]
        time <- 1:end[ii]
        lam <- calc.lambda(time, b0=b0, b1=b1)
        cbind( time=max(time) - rev(time), lam, group=ii, sim=kk)
    }) # end coef pair
    res <- data.frame( do.call(rbind, cPair) )
    res
}) %>% rbind.fill # end group


#
#       END \\ SUMMARIZE MODEL RESULTS
#


#
#       PLOT MODEL FIT
#


# relabel index to original scale
both$year <- both$index + min(dataRAW$year)

# set up plotting data for group panels
obs <- filter(both, sim==0) # subset to observed series
sims <- filter(both, sim!=0) %>% # subset a sample of simmed series
    split( . , .$group) %>% # group by group
    lapply( . , function(oo){ # for each group
        ids <- sample(unique(oo$sim), 200) # sample 200sims
        oo[ oo$sim %in% ids, ] # subset sims
    } ) %>% rbind.fill # bind together

# set up facet labels
labels <- as.character(mapping$groupname)
names(labels) <- mapping$group

# plot cumulative counts
P <- ggplot( ) +
    geom_path(data=sims, aes(x=year, y=cml_value, group=sim), 
        col="skyblue2", alpha=0.1) +
    geom_path(data=obs, aes(x=year, y=cml_value)) +
    facet_wrap(~group, scales="free_y", labeller=as_labeller(labels) ) +
    ylab("Number of Species") + 
    xlab("Year of Description")
ggsave(P, file="output/cumulative_fit.pdf", width=10, height=6)

# plot per year counts
P <- ggplot( ) +
    geom_path(data=sims, aes(x=year, y=value, group=sim), 
        col="skyblue2", alpha=0.1) +
    geom_path(data=obs, aes(x=year, y=value)) +
    facet_wrap(~group, scales="free_y", labeller=as_labeller(labels) ) +
    ylab("Number of Species") + 
    xlab("Year of Description")
ggsave(P, file="output/count_fit.pdf", width=10, height=6)


# set up plotting data for long-term trends
lambda$year <- lambda$time + min(dataRAW$year)
sims <- filter(lambda, sim!=0) %>% # subset a sample of simmed series
    split( . , .$group) %>% # by group
    lapply( . , function(oo){ # for each group
        ids <- sample(unique(oo$sim), 200) # sample 200sims
        oo[oo$sim %in% ids, ] # subset sims
    } ) %>% rbind.fill # bind
# mean line per group
mu_sim <- sims %>%
    group_by(group, year) %>%
    dplyr::summarize(mu=mean(lam)) %>%
    data.frame()

# Plot long term trends
P <- ggplot( ) +
    geom_path(data=obs, aes(x=year, y=value)) +
    geom_line(data=sims, aes(x=year, y=lam, group=sim), col="skyblue2", 
        alpha=0.1) +
    geom_line(data=mu_sim, aes(x=year, y=mu), col="royalblue") +
    facet_wrap(~group, scales="free_y", labeller=as_labeller(labels) ) +
    ylab("Number of Species") + 
    xlab("Year of Description")
ggsave(P, file="output/regression.pdf", width=10, height=6)


#
#       END \\ PLOT MODEL FIT
#



#
#       SUMMARIZE FORECAST
#

# cumulative series for sampled data
forsim <- lapply( seq(length(forecast)), function(jj) { # each sim
    lapply( seq(length(data$end)), function(ii) { # each prov
        data.frame(
            index=((ncol(data$counts) + 1) : (ncol(data$counts)+
                length(forecast[[1]][[1]]))) + min(dataRAW$year) - 1, 
            value=cumsum(forecast[[jj]][[ii]]), group=ii, sim=jj
        )
    } ) %>% rbind.fill
}) %>% rbind.fill

# for each prov
# cumsum across each sim
# find the mean cumsum
# find the CI of cumsums

fore.table <- split(forsim, forsim$group) %>% # by group
    lapply( . , function(gg) { # for each group
        group <- unique(gg$group)
        vv <- lapply( split(gg, gg$sim), function(oo) { # for each sim
            cs <- max(oo$value)
            oo[ nrow(oo), ]$value <- cs
            oo[ nrow(oo), ]
        }) %>% rbind.fill
        fore.mu <- round( mean(vv$value), 0)  + # summarize mean expected
            max(obs[ obs$group == group, "value"]) # add mean expected to currents counts
        fore.lower <- round( quantile(vv$value, 0.1), 0)  + # summarize lower expected
            max(obs[ obs$group == group, "value"]) # add mean expected to currents counts    
        fore.upper <- round( quantile(vv$value, 0.9), 0)  + # summarize upper expected
            max(obs[ obs$group == group, "value"]) # add mean expected to currents counts
        data.frame(group=group, fore.mu, fore.lower, fore.upper)    
    }) %>% rbind.fill

# merge to Results table
RESULTS <- merge(results, fore.table, by="group") %>%
    arrange(desc(observed_species))
write.csv(RESULTS, file="output/results.csv", row.names=FALSE)

#
#       END \\ SUMMARIZE FORECAST
#