# run all of the prep- files in the right order

# adjust options.R, if needed

# if prep dir is present, then an old run probably exists that
# should be removed first
stopifnot(!file.exists(getOption('dropviz.prep.dir')))

files <- c('global','tSNE','metacells','components','marker','expr')
lapply(files, function(fn) {
  cat("Sourcing ", fn, "\n")
  source(paste0("prep-",fn,".R"))
})

