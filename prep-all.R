# run all of the prep- files in the right order

options(dropviz.experiments='exp_sets.full.txt',
        dropviz.prep.dir='staged', # where to write prep'd data. Move this to the live location when done
        dropviz.prep.cache=FALSE)  # use existing files. Usually set this to FALSE during a prep to force since there are interdependencies that affect files such as the factor values

# if prep dir is present, then an old run probably exists that
# should be removed first
stopifnot(!file.exists(getOption('dropviz.prep.dir')))

files <- c('tSNE','metacells','components','marker','expr')
lapply(files, function(fn) {
  cat("Sourcing ", fn, "\n")
  source(paste0("prep-",fn,".R"))
})

