# run all of the prep- files in the right order

# if cache dir is present, then an old run probably exists that
# should be removed first
stopifnot(!file.exists("cache"))

files <- c('tSNE','expr','metacells','components','marker')
lapply(files, function(fn) {
  cat("Sourcing ", fn, "\n")
  source(paste0("prep-",fn,".R"))
})

