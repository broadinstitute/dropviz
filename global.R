# library(strict) - issues spurious warnings that I can't track down
library(shiny)
library(shinyjs)
library(shinyjqui)
library(plyr)
library(tidyr)
library(dplyr)
library(glue)
library(ggthemes)
library(digest)
library(ggplot2)

source("shared.R")

log.reactive <- function(...) {
  if (getOption("log.reactive", default=FALSE)) write.log(...)
}

# gracefully fails returning a NULL if not in a reactive environment
shiny.progress <- function(message=NULL) {
  tryCatch({
    p <- shiny::Progress$new()
    p$set(message=message)
    p
  }, error=function(err) NULL)
}

write.func.body <- function(fn, file) {
  fn.lines <- capture.output(print(fn))
  # skip anonymous function def (first line), environment label (last line) and closing brace (penultimate line)
  writeLines(fn.lines[2:(length(fn.lines)-2)], file)
}

send.zip <- function(fn, fname, zipfile) {
  require(utils)
  zip.dir <- tempdir()
  fn.env <- environment(fn)
  vars <- ls(environment(fn))
  cwd <- getwd()
  setwd(zip.dir)
  zip.files <- paste0(fname, c(".Rdata", ".R"))
  attach(fn.env)
  save(list=vars, file=zip.files[1])
  detach()
  write.func.body(fn, file=zip.files[2])
  zip(zipfile, zip.files)
  setwd(cwd)
}

if (file.exists(glue("{prep.dir}/markers/top_genes.RDS"))) {
  top.genes <<- sort(readRDS(glue("{prep.dir}/markers/top_genes.RDS")))
}

## read gene symbols and descriptions
if (file.exists(glue("{prep.dir}/markers/gene.dict.RDS"))) {
  gene.dict <<- readRDS(glue("{prep.dir}/markers/gene.dict.RDS"))
}

# returns empty string for failed matches
gene.lookup <- function(name, pos) {
  if (is.null(gene.dict)) {
    gene.dict <<- mk.gene.dict()
  }
  if (is.null(name) || length(name)==0) {
    character(0)
  } else {
    lapply(tolower(name), function(x) if (exists(x,envir=gene.dict)) gene.dict[[x]][pos] else "") %>% unlist
  }
}

gene.symbol <- function(name) {
  gene.lookup(name, 1)
}

gene.desc <- function(name) {
  gene.lookup(name, 2)
}

stopifnot(file.exists(glue("{prep.dir}/globals.Rdata")))
load(glue("{prep.dir}/globals.Rdata"))
