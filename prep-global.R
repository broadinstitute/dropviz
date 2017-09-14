# save project-wide data to cache/globals.Rdata:
#
# experiments - read from exp_sets.txt or option dropviz.experiments
# cell.types - cleaned, concatenation of all cluster_sheets/ CSVs.
# cluster.names_, subcluster.names_ - meaningful names for (sub)clusters. usually use reactives instead of these directly
# components - cleaned, concatenation of all curation_sheets/ CSVs

source("shared.R")

library(plyr)
library(dplyr)
library(glue)

options(stringsAsFactors=FALSE) # plyr methods return data.frame

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
        cluster_name <- glue("{paste0(unique(df$class),collapse='/')}.{first.class_marker}")
      } else {
        if (all(df$class==first.class)) {
          cluster_name <- first(df$class)
        } else {
          # oh man, these are ugly clusters. Example is frontal cortex 11:
          # A tibble: 4 x 8
          # exp.label subcluster cluster         region      class class_marker                  full_name
          # <fctr>     <fctr>  <fctr>          <chr>      <chr>        <chr>                      <chr>
          # 1 F_GRCm38.81.P60Cortex_noRep5_FRONTALonly       11-1      11 Frontal Cortex  Microglia         C1qb     Microglia.C1qb.Tmem119
          # 2 F_GRCm38.81.P60Cortex_noRep5_FRONTALonly       11-2      11 Frontal Cortex     Neuron      Slc17a6        Neuron.Slc17a6.Reln
          # 3 F_GRCm38.81.P60Cortex_noRep5_FRONTALonly       11-3      11 Frontal Cortex Macrophage         C1qb       Macrophage.C1qb.Mrc1
          # 4 F_GRCm38.81.P60Cortex_noRep5_FRONTALonly       11-4      11 Frontal Cortex  Microglia         C1qb Microglia.C1qb.Tmem119-Fos
          #
          # Use class.marker
          cluster_name <- paste(unlist(dlply(df, .(class_marker), function(df2) paste0(paste(unique(df2$class),collapse='/'),'.',first(df2$class_marker)))),collapse='+')
          write.log("Ugly: ",cluster_name)
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
experiments.fn <- getOption('dropviz.experiments', default='exp_sets.txt')
experiments <- as_tibble(read.delim(experiments.fn, header = TRUE, stringsAsFactors = FALSE)) %>%
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
    write.log(cluster.sheets.dir,'/',fn)
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


# Create "biologically meaningful" names for cluster and subcluster
cluster.names_ <- mk.cluster.names(cell.types_)
subcluster.names_ <- mk.subcluster.names(cell.types_)
cell.types <- select(cell.types_, -region, -full_name, -common_name, -class_marker) 


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

save(experiments, cell.types, cluster.names_, subcluster.names_, components, file=glue("{prep.dir}/globals.Rdata"))
