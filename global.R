# library(strict) - issues spurious warnings that I can't track down
library(shiny)
library(shinyjs)
library(shinyjqui)
library(plyr)
library(dplyr)
library(glue)

options(stringsAsFactors=FALSE) # plyr methods return data.frame
options(bitmapType = 'cairo') # https://stackoverflow.com/a/17955000/86228

# write either warning or message to console with a time stamp
write.log <- function(..., warn=FALSE) {
  level.func <- ifelse(warn, warning, message)
  #  level.func(Sys.time(),": ",glue(.envir=parent.frame(2), ...))  # FIXME
  level.func(Sys.time(),": ", ...)
}

log.reactive <- function(...) {
  if (getOption("log.reactive", default=FALSE)) write.log(...)
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

readRDS <- function(file, ...) {
  # for performance tracking, log when going to disk
  write.log("Reading ",file)
  base::readRDS(file, ...)
}



MAX_REGIONS <- 4

prep.dir <- 'cache'               # destination for all prepared data
cache.dir <- 'www/cache'  # destination for renderCacheImage()
suppressWarnings(dir.create(prep.dir, recursive = TRUE))
suppressWarnings(dir.create(cache.dir, recursive = TRUE))

# Replaces strings of form N-Month with M-N.
# This is a hack to fix Excel converting M-N into a date format.
patch.cluster.number <- function(cn) {
  month.map <- c(Jan=1,Feb=2,Mar=3,Apr=4,May=5,Jun=6,Jul=7,Aug=8,Sep=9,Oct=10,Nov=11,Dec=12)
  month.names <- names(month.map)
  invisible(lapply(month.names, function(m) {
    month.regex <- paste0('([0-9]+)-',m)
    month.repl <- paste0(month.map[[m]],"-\\1")
    cn <<- sub(month.regex, month.repl, cn)  
  }))
  cn
}

# TODO: I later discovered that curation_sheets/*.cluster_class contains a cluster_name for each cluster.
# Use that instead!
#
# In the cell.types tibble, each class cluster can be assigned a name by choosing the class and marker
# E.g.
# 11-1 Astrocyte.Gja1.Myoc
# 11-2 Astrocyte.Gja1.Vegfa
# yields Astrocyte.Gja1
# Sometimes there is a common marker, but the class name differs:
# Macrophage.C1qb.Mrc1
# Microglia.C1qb.Tmem119
# In those cases, choose the class names like "Macrophage/Microglia:C1qb"
mk.cluster.names <- function(ct) {
  
  ddply(ct, .(exp.label, cluster), function(df) {
    first.class <- first(df$class)
    first.class_marker <- first(df$class_marker)
    
    # The usual case
    if (all(df$class==first.class) && all(df$class_marker==first.class_marker)) {
      cluster_name <- glue("{first.class}.{first.class_marker}")
    } else {
      if (all(df$class_marker==first.class_marker)) {
        # The Macrophage/Microglia case
        cluster_name <- glue("{paste0(df$class,collapse='/')}.{first.class_marker}")
      } else {
        if (all(df$class==first.class)) {
          cluster_name <- first(df$class)
        } else {
          write.log("Stuck in mk.cluster.names on ", first(df$exp.label),' ',first(df$cluster))
        }
      }
    }
    data.frame("cluster_name"=cluster_name)
  }) %>% as_tibble
}

# Make a subcluster_name that is the full_name unless there is a common_name. If so, use that
mk.subcluster.names <- function(ct) {
  sc.names <- with(ct, tibble(exp.label=exp.label, subcluster=subcluster,
                                  full_name=full_name,
                                  subcluster_name=ifelse(is.na(common_name)|common_name=="", full_name, common_name)))
}


# read a config file that lists the names and locations of each experiment set, e.g.
#
# exp.label exp.title                                           exp.dir
# 1 GRCm38.81.P60Striatum  Striatum /cygwin64/home/dkulp/data/F_GRCm38.81.P60Striatum
# 2 GRCm38.81.P60Thalamus  Thalamus /cygwin64/home/dkulp/data/F_GRCm38.81.P60Thalamus
experiments <- as_tibble(read.delim("exp_sets.txt", header = TRUE, stringsAsFactors = FALSE)) %>%
  mutate(exp.label=as.factor(exp.label))
experiments$base <- basename(experiments$exp.dir)


# read the cluster_sheets/cluster_N.csv files and concatenate them all together.
#
# A tibble: 6 × 13
# exp.label                  tissue cluster_number  class class_marker type_marker                  full_name
# <fctr>                   <chr>          <chr>  <chr>        <chr>       <chr>                      <chr>
# 1 GRCm38.81.P60Striatum GRCm38.81.P60Cerebellum          1-1 Neuron      Slc17a7      Gabra6      Neuron.Slc17a7.Gabra6
# 2 GRCm38.81.P60Striatum GRCm38.81.P60Cerebellum          2-2 Neuron      Slc17a7  Gabra6-Fos  Neuron.Slc17a7.Gabra6-Fos
# 3 GRCm38.81.P60Striatum   GRCm38.81.P60Striatum          10-1 Neuron     Gad1Gad2    Drd1-Fos   Neuron.Gad1Gad2.Drd1-Fos
# 4 GRCm38.81.P60Striatum   GRCm38.81.P60Striatum          10-2 Neuron     Gad1Gad2  Drd1-Lypd1 Neuron.Gad1Gad2.Drd1-Lypd1
#
cell.types_ <- as_tibble(ddply(experiments, .(exp.label), function(exp) {
  cluster.sheets.dir <- sprintf("%s/cluster_sheets",exp$exp.dir)
  ldply(list.files(cluster.sheets.dir, "Cluster_[0-9]+.csv"), function(fn) {
    # FIXME: warning suppressed because last line is missing CR
    csv <- suppressWarnings(read.csv(sprintf("%s/%s", cluster.sheets.dir, fn), fill=TRUE, stringsAsFactors=FALSE)) 
    
    # Add the experiment info and remove some unused columns
    csv <- dplyr::mutate(csv, subcluster=patch.cluster.number(cluster_number),
                         exp.label=exp$exp.label,
                         exp.title=exp$exp.title,
                         cluster=sub("([0-9]+)-[0-9]+", "\\1", subcluster)) %>% 
      select(exp.label, subcluster, cluster, region=exp.title, class, class_marker, full_name, common_name) %>%
      mutate(cluster=as.factor(cluster), subcluster=as.factor(subcluster))
  })
})) %>% mutate(exp.label=as.factor(exp.label))

#cell.types$cluster_disp <- sprintf("%s %s (%s)", cell.types$class, cell.types$cluster, cell.types$exp.title)


# Create "biologically meaningful" names for cluster and subcluster
cluster.names_ <- mk.cluster.names(cell.types_)
subcluster.names_ <- mk.subcluster.names(cell.types_)
cell.types <- select(cell.types_, -region, -full_name, -common_name, -class_marker) 

# # returns a tibble of class and marker gene
# class.markers <- {
#   markers.split <- strsplit(cell.types$class_marker,'-')
#   ldply(1:length(markers.split), function(n) {
#     tibble(class_marker=markers.split[[n]], 
#            class_markers=cell.types$class_marker[n], class=cell.types$class[n])
#   }, .id=NULL) %>% unique
# }


# type.markers <- {
#   markers.split <- strsplit(cell.types$type_marker,'-')
#   ldply(1:length(markers.split), function(n) {
#     tibble(type_marker=markers.split[[n]],
#            type_markers=cell.types$type_marker[n], full_name=cell.types$full_name[n])
#   }, .id=NULL) %>% unique
# }

# all genes that are differentially expressed according to the most liberal criteria
# in all experiments and all comparison methods
# there are too many (16,700) gene names to load them all in the client
# for testing, just sample a few or use the annotated markers
#all.genes <- sample(readRDS("markers/all.genes.RDS"), 1000)
#all.genes <- unique(unlist(strsplit(type.markers$type_marker,'\\.')))
markers.fn <- glue("{prep.dir}/markers/pval-200.genes.RDS")
if (file.exists(markers.fn)) {
  all.genes <- sort(readRDS(markers.fn))
}

components <- as_tibble(ddply(experiments, .(exp.label), function(exp) {
  curation.sheets.dir <- sprintf("%s/curation_sheets",exp$exp.dir)
  ldply(list.files(curation.sheets.dir, "AnnotationSheet_Subcluster_[0-9]+.*csv"), function(fn) {
    cluster <- regmatches(fn, regexpr("[0-9]+",fn))
    csv <- suppressWarnings(read.csv(sprintf("%s/%s", curation.sheets.dir, fn), fill=TRUE, stringsAsFactors=FALSE)) %>%
      mutate(status=factor(status),
             use_for_clustering=ifelse(use_for_clustering=='Y',TRUE,ifelse(use_for_clustering=='N',FALSE,NA)))
    cbind(data.frame(cluster=factor(cluster,levels=levels(cluster.names_$cluster))), csv)
  })
}))

# read gene descriptions
gene.desc.dict <- NULL
gene.desc <- function(name) {
  if (is.null(gene.desc.dict)) {
    gene.descriptions <- read.delim("data/gene_descriptions.txt.gz") %>% mutate(Description=sub(' \\[.*','',Description))
    gene.desc.list <- setNames(lapply(1:nrow(gene.descriptions), function(i) gene.descriptions$Description[i]), gene.descriptions[['Associated.Gene.Name']])
    gene.desc.dict <<- list2env(gene.desc.list)
  }
  if (is.null(name) || length(name)==0) {
    character(0)
  } else {
    lapply(lapply(name, function(x) gene.desc.dict[[x]]), function(n) if (is.null(n)) "" else n) %>% unlist
  }
}
