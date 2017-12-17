source("options.R")

options(stringsAsFactors=FALSE) # plyr methods return data.frame
options(bitmapType = 'cairo') # https://stackoverflow.com/a/17955000/86228

prep.dir <- getOption('dropviz.prep.dir', default='staged')               # destination for all prepared data
cache.dir <- 'www/cache'  # destination for renderCacheImage()
suppressWarnings(dir.create(prep.dir, recursive = TRUE))
stopifnot(file.exists(prep.dir))
suppressWarnings(dir.create(cache.dir, recursive = TRUE))
stopifnot(file.exists(cache.dir))

# write either warning or message to console with a time stamp
write.log <- function(..., warn=FALSE) {
  level.func <- ifelse(warn, warning, message)
  #  level.func(Sys.time(),": ",glue(.envir=parent.frame(2), ...))  # FIXME
  level.func(Sys.time(),": ", ...)
}

readRDS <- function(file, ...) {
  # for performance tracking, log when going to disk
  write.log("Reading ",file)
  base::readRDS(file, ...)
}

saveRDS <- function(object, file, ...) {
  write.log("Writing ", substitute(object), " to ", file)
  base::saveRDS(object, file, ...)
}

