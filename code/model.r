# fit model to data
library(rstan)

# initial data
files <- dir('data/dump', pattern = 'count_info.data.R',
             full.names = TRUE)
data <- read_rdump(files)

# fit model
fit <- stan( file="code/zip_count.stan",
                   data=data, 
                   chains=4, 
                   warmup=500,
                   iter=1000,
                   init=0,
                   thin=5,
                   cores=4,
                   verbose=TRUE)
save(fit, file="data/dump/fit.data")
