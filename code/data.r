# build data object for model fit
args <- commandArgs(trailingOnly=TRUE)
OFFSET <- as.logical(args[1])

library(plyr)
library(rstan)

# change as necessary
data <- read.csv("data/data.csv")

duration <- seq(from = min(data$year), to = max(data$year))
group <- split(data, data$group)

counts <- llply(group, function(x) table(x$year))
counts <- llply(counts, function(x) {
                pos <- rep(0, length(duration))
                names(pos) <- duration
                pos[duration %in% names(x)] <- x
                pos})
count.matrix <- Reduce(cbind, counts)  # named per year in a group


# actually fit a model
# got to get data in stan format
nyear <- nrow(count.matrix)
jgroup <- ncol(count.matrix)
npred <- 1

starts <- apply(count.matrix, 2, function(x) min(which(x != 0)))
count.matrix[starts[1], 1]
data <- list()
for(ii in seq(jgroup)) {
  long <- length(count.matrix[, ii])
  data[[ii]] <- list(counts = count.matrix[seq(from = starts[ii], 
                                               to = long), ii],
                     N = length(seq(from = starts[ii], to = long)))
}
years.named <- llply(data, function(x) names(x$counts))


# switch it up, making first col the first year things named
cc <- list()
len <- laply(data, function(x) length(x$counts))
for(ii in seq(length(data))) {
  cc[[ii]] <- c(data[[ii]]$counts, rep(0, (max(len) - len[ii])))
}
cc <- Reduce(rbind, cc)
colnames(cc) <- NULL
rownames(cc) <- NULL
N <- ncol(cc)
P <- nrow(cc)

# number of publications per year
npub <- llply(group, function(x) 
              table(unique(x[, c('species_authority', 'year')])$year))
pub.matrix <- array(0, dim = dim(count.matrix))
rownames(pub.matrix) <- rownames(count.matrix)
for(ii in seq(length(npub))) {
  pub.matrix[rownames(pub.matrix) %in% names(npub[[ii]]), ii] <- npub[[ii]]
}


if( ! OFFSET ){  # set offset to zero
    pub.matrix[] <- 0
}

data <- list(N = N, P = P, end = len, counts = cc, off = t(pub.matrix))

with( data, {stan_rdump(list = c('N', 'P', 'end', 'counts', 'off'),
    file = 'data/dump/count_info.data.R')} )