# visualize results from model fit and forecasting
args <- commandArgs(trailingOnly=TRUE)
OFFSET <- as.logical(args[1])
# OFFSET <- TRUE

library(dplyr)
library(plyr)
library(rstan)
library(ggplot2)


#    LOAD DATA --------------------------------------------------------------------------------------------

dataRAW <- read.csv("data/data.csv") # original data
data <- read_rdump("data/dump/count_info.data.R") # initial model data


load("data/dump/fit.data") # zero inflated fits --loads as fit
zips <- fit # reassign to zips
rm(fit) # remove from memory
load("data/dump/post.data") # posterior samples -- loads as allsim
load("data/dump/forecast.data") # forecast data -- loads as forecast

# map model indeces to original variables
mapping <- data.frame( 
    groupname=as.character(unique(dataRAW$group)),
    group=as.numeric(unique(dataRAW$group)))

#   END \\ LOAD --------------------------------------------------------------------------------------------


#   PREPARE CUMULATIVE SERIES -------------------------------------------------------------------

# cumulative series for observed data
cumm <- lapply( seq(data$P), function(ii) { # each prov
    data.frame(
        index=1:data$end[ii],
        value=data$counts[ii, ],
        cml_value=cumsum(data$counts[ii, ]),
        off=data$off[ii, ],
        group=ii
    )
}) %>% rbind.fill
cumm$sim <- 0 # set sim to 0 to mark observed data

# cumulative series for simmed data
cummsim <- lapply( seq(length(allsim)), function(jj) { # each sim
    lapply( seq(length(data$end)), function(ii) { # each prov
        data.frame(
            index=1:data$end[ii], 
            value=allsim[[jj]][[ii]],
            cml_value=cumsum(allsim[[jj]][[ii]]),
            off=data$off[ii, ],
            group=ii,
            sim=jj
        )
    } ) %>% rbind.fill
}) %>% rbind.fill

Z <- rbind(cumm, cummsim) # combine the observed and simmed series
Z$year <- Z$index + min(dataRAW$year) - 1 # add original year back

#   END \\ PREPARE CUMULATIVE SERIES ---------------------------------------------------------


#   SUMMARIZE MODEL RESULTS --------------------------------------------------------------------

sumy <- filter(Z, sim !=0 & index==max(index)) %>%
  group_by(group) %>%
  dplyr::summarize(
    med=median(cml_value),
    lower=quantile(cml_value, probs=0.1),
    upper=quantile(cml_value, probs=0.9)
    )
obs.count <- filter(Z, sim==0 & index==max(index)) %>%
  group_by(group) %>%
  dplyr::summarize(count=cml_value)

# outs$coef[2] is the coefficient that estimates the long-term trend in
# description rate, it's in log units so need exponentiate it
# where 0 means no trend, positive means increasing description
#  through time and negative is decreasing...
# `coef` is an array where 
# [posterior sample [1:4000], province number[1:16], coefficient[1:2] ]
coef <- rstan::extract(zips, par="coef")[[1]]
# pull out the second coefficient for each province across each 
# posterior sample
cf2 <- apply(coef, 1, function(i) i[ , 2] )
mean_cf2 <- data.frame(provid=1:nrow(cf2),
    slowdown=apply(cf2, 1, mean) )
CI80_cf2 <- data.frame(provid=1:nrow(cf2), 
    slowdown=t( apply(cf2, 1, quantile, probs=c(0.1, 0.9))) )

# compile results into table
results <- data.frame(
    group=mapping[ match(as.numeric(rownames(obs.count)), mapping[, 2]), 2 ] , 
    groupname=mapping[ match(as.numeric(rownames(obs.count)), mapping[, 2]), 1 ] , 
    observed_species=obs.count$count, 
    expected_median=round(sumy$med, 0), 
    expected_CI_lower=round(sumy$lower, 0), 
    expected_CI_higher=round(sumy$upper, 0),
    slowdown=round( signif(mean_cf2$slowdown, 3 ), 4 ),
    slowdown_CI_lower=round( signif(CI80_cf2[,2], 3), 4),
    slowdown_CI_higher=round( signif(CI80_cf2[,3], 3), 4)
,row.names=NULL) %>% arrange( slowdown )


#   Long-term regression lines
calc.lambda <- function(x, b0, b1) {
  exp(b0 + b1 * x)
}

# pull out the second coefficient for each province across each posterior sample
prov.cf1 <- apply(coef, 1, function(i) i[ ,1])
prov.cf2 <- apply(coef, 1, function(i) i[ ,2])

lambda <- lapply( seq(data$P), function(ii) { # for each group
    cf1 <- prov.cf1[ii, ]
    cf2 <- prov.cf2[ii, ]
    time <- 1:data$end[ii]
    cPair <- lapply(1:length(cf1), function(kk){ # for each coefficient pair
        b0 <- cf1[kk]
        b1 <- cf2[kk]
        lam <- calc.lambda(time, b0=b0, b1=b1)
        cbind( time=max(time) - rev(time), lam, group=ii, sim=kk)
    }) # end coef pair
    res <- data.frame( do.call(rbind, cPair) )
    res
}) %>% rbind.fill # end group
lambda$year <- lambda$time + min(dataRAW$year) # add original year back

#       END \\ SUMMARIZE MODEL RESULTS ------------------------------------------------------


#       PLOT MODEL FIT -----------------------------------------------------------------------------------

# set up plotting data for group panels
obs <- filter(Z, sim==0) # subset to observed series
sims <- filter(Z, sim!=0) %>% # subset a sample of simmed series
    split( . , .$group) %>% # group by group
    lapply( . , function(oo){ # for each group
        ids <- sample(unique(oo$sim), 200) # sample 200sims
        oo[ oo$sim %in% ids, ] # subset sims
    } ) %>% rbind.fill # bind together
# only keep sims less that 4 times the max observed value
goodsims <- filter(sims, year==max(year) & cml_value < max(obs$cml_value)*4) %>%
    select(group, sim)
sims <- inner_join(sims, goodsims, by=c("group", "sim"))

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
    geom_path(data=obs, aes(x=year, y=(value/(off+1)))) +
    geom_line(data=sims, aes(x=year, y=lam, group=sim), col="skyblue2", 
        alpha=0.1) +
    geom_line(data=mu_sim, aes(x=year, y=mu), col="royalblue") +
    facet_wrap(~group, scales="free_y", labeller=as_labeller(labels) ) +
    ylab("Number of Species") + 
    xlab("Year of Description")
ggsave(P, file="output/regression.pdf", width=10, height=6)


#       END \\ PLOT MODEL FIT -------------------------------------------------------------------------


#       SUMMARIZE FORECAST --------------------------------------------------------------------------

# cumulative series for sampled data
forsim <- lapply( seq(length(forecast)), function(jj) { # each sim
    lapply( seq(data$P), function(ii) { # each prov
        data.frame(
            index=((ncol(data$counts) + 1) : (ncol(data$counts) +
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
        # gg <- split(forsim, forsim$group)[[1]]
        group <- unique(gg$group)
        vv <- lapply( split(gg, gg$sim), function(oo) { # for each sim
            cs <- max(oo$value)
            oo[ nrow(oo), ]$value <- cs
            oo[ nrow(oo), ]
        }) %>% rbind.fill
        fore.mu <- round( mean(vv$value), 0)  + # summarize mean expected
            max(obs[ obs$group == group, "cml_value"]) # add mean expected to currents counts
        fore.lower <- round( quantile(vv$value, 0.1), 0)  + # summarize lower expected
            max(obs[ obs$group == group, "cml_value"]) # add mean expected to currents counts    
        fore.upper <- round( quantile(vv$value, 0.9), 0)  + # summarize upper expected
            max(obs[ obs$group == group, "cml_value"]) # add mean expected to currents counts
        data.frame(group=group, fore.mu, fore.lower, fore.upper)    
    }) %>% rbind.fill

# merge to Results table
RESULTS <- merge(results, fore.table, by="group") %>%
    arrange(desc(observed_species))
write.csv(RESULTS, file="output/results.csv", row.names=FALSE)

#       END \\ SUMMARIZE FORECAST ---------------------------------------------------------------