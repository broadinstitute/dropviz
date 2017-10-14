# run all of the prep- files in the right order

# adjust options.R, if needed


files <- c('tSNE','metacells','components','marker','expr')
lapply(files, function(fn) {
  cat("Sourcing ", fn, "\n")
  source(paste0("prep-",fn,".R"))
})

