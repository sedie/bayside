# wrap for running bayside model
args <- commandArgs(trailingOnly=TRUE)
OFFSET <- args[1]  # TRUE = include publications in model, FALSE = do not include publications
ftime <- args[2]  # number of years into the future to forecast species richness (integer; 0 -> INF)

# run data
system(paste0("R CMD BATCH --vanilla '--args ", OFFSET, "' code/data.r"))
# run model
system("R CMD BATCH --vanilla code/model.r")
# run post
system("R CMD BATCH --vanilla code/post.r")
# run forecast
system(paste0("R CMD BATCH --vanilla '--args ", ftime, "' code/forecast.r "))
# run plot
system(paste0("R CMD BATCH --vanilla '--args ", OFFSET, "' code/plot.r"))