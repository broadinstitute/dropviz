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

stack.trace <- capture.output(traceback(0,1))
if (any(grepl("serviceApp|runApp", stack.trace, perl=TRUE))) {
  # hack. Checking the call stack for runApp is the best way I can find to determine if
  # I'm running as a server or not.
  message("Removing local reactive overrides, if any")
  suppressWarnings(try(rm(reactive, reactiveValues, observeEvent, output, onRestore), silent=TRUE))
} else {
  message("Running in interactive, non-Shiny environment")
  reactive <- function(x, env=parent.frame(), ...) exprToFunction(x, env=env)
  reactiveValues <- list
  output <- list()
  onRestore <- function(...) {}
  observeEvent <- function(...) {}
}

# loads the bookmark into the input variable in the caller's environment
load.mark <- function(id) assign('input',readRDS(glue("shiny_bookmarks/{id}/input.rds")), 1)

source("shared.R")

# gracefully fails returning a NULL if not in a reactive environment
shiny.progress <- function(message=NULL) {
  tryCatch({
    p <- shiny::Progress$new()
    p$set(message=message)
    p
  }, error=function(err) NULL)
}

# take the environment in the function fn and save to fname.Rdata
# this function used to also print the function body, but it appears inconsistent across platforms.
# add others after removing first and last line.
# zip them all up.
send.zip <- function(fn, fname, zipfile, others=character(0)) {
  require(utils)
  zip.dir <- tempdir()
  fn.env <- environment(fn)
  vars <- ls(environment(fn))
  file.copy(others, zip.dir)
  cwd <- getwd()
  setwd(zip.dir)
  zip.files <- c(paste0(fname, ".Rdata"), basename(others))
  attach(fn.env)
  save(list=vars, file=zip.files[1])
  detach()
  zip(zipfile, zip.files)
  setwd(cwd)
}

if (file.exists(glue("{prep.dir}/markers/top_genes.RDS"))) {
  top.genes <<- sort(readRDS(glue("{prep.dir}/markers/top_genes.RDS")))
  if (!file.exists(file.path('www/top.genes.json'))) {
    writeLines(paste('{ "genes": [\n',
                     paste0('{ "value": "', top.genes, '", "label": "', top.genes, '" }', collapse=",\n"),
                     '\n]}'), "www/top.genes.json")
  }
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

# #34 replace old names
# Endothelial_Tip -> Fibroblast-Like
# Endothelial_Stalk -> Endothelial
# Ependyma.+ -> Ependyma

cluster.names_ <- 
  mutate(cluster.names_, 
                         class=sub('Ependyma.+', 'Ependyma', sub('Endothelial_Stalk', 'Endothelial', sub('Endothelial_Tip', 'Fibroblast-Like', class))),
                         cluster_name=sub('Ependyma.+', 'Ependyma', sub('Endothelial_Stalk', 'Endothelial', sub('Endothelial_Tip', 'Fibroblast-Like', cluster_name)))
)

subcluster.names_ <- mutate(subcluster.names_,
                            full_name=sub('Ependyma.+', 'Ependyma', sub('Endothelial_Stalk', 'Endothelial', sub('Endothelial_Tip', 'Fibroblast-Like', full_name))),
                            subcluster_name=sub('Ependyma.+', 'Ependyma', sub('Endothelial_Stalk', 'Endothelial', sub('Endothelial_Tip', 'Fibroblast-Like', subcluster_name))))

components <- mutate(components,
                     cell_class=gsub('Ependyma.+', 'Ependyma', gsub('Endothelial_[Ss]talk', 'Endothelial', gsub('Endothelial_[Tt]ip', 'Fibroblast-Like', cell_class))))

# #36 replace exp.title
experiments <- mutate(experiments, exp.title=sub('Ento Peduncular', 'Entopeduncular', exp.title))
