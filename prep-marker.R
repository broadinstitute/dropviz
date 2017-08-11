# combine the gene expression results for all clusters and sub-cluster with respect to global and cluster
# write the results in classMarkers.RDS typeMarkers.RDS

library(plyr)
library(dplyr)
library(glue)

source("global.R")
markers.dir <- glue("{prep.dir}/markers")

read.clusterMarkers <- function(exp, kind='cluster') {
    if (kind=='subcluster') {
        fn <- '_clusterMarkers.RDS'
    } else {
        fn <- "_FirstRound_clusterMarkers.RDS"
    }
    clusterMarkers <- readRDS(sprintf("%s/%s%s",exp$exp.dir,exp$base,fn))

    as_tibble(mutate(ldply(clusterMarkers, function(df) {
      df
    }, .id="cluster"), cluster=as.character(cluster)))
}


# # A tibble: 42,004 x 7
# cluster    GENE FOLDch  pval pct.1 pct.2             exp.label
# <chr>   <chr>  <dbl> <dbl> <dbl> <dbl>                 <chr>
# 1       1 Rarres2 51.687     0  0.98  0.02 GRCm38.81.P60Striatum
# 2       1 Ccdc153 45.053     0  0.94  0.00 GRCm38.81.P60Striatum
# 3       1     Dbi 31.138     0  0.99  0.36 GRCm38.81.P60Striatum
# 4       1    Nnat 30.840     0  0.97  0.12 GRCm38.81.P60Striatum
readall.clusterMarkers <- function(kind) {
  as_tibble(ddply(experiments, .(exp.label), function(exp) {
    mutate(read.clusterMarkers(exp, kind), exp.label=exp$exp.label)
  }))
}

clusterMarkers <- readall.clusterMarkers('cluster') %>% mutate(cluster=factor(cluster,levels=levels(cluster.names_$cluster)))
subclusterMarkers <- readall.clusterMarkers('subcluster') %>% mutate(cluster=factor(cluster,levels=levels(subcluster.names_$subcluster)))
dir.create(markers.dir, recursive = TRUE)
saveRDS(clusterMarkers, file=glue("{markers.dir}/clusterMarkers.RDS"))
saveRDS(subclusterMarkers, file=glue("{markers.dir}/subclusterMarkers.RDS"))
all.genes <- unique(c(clusterMarkers$GENE, subclusterMarkers$GENE))

########################################
## For the pairwise, I store each RDS separately
## per experiment and cluster because the user can only choose
## to compare to another class or type within an experiment
## and individual access is faster.
read.pairwiseMarkers <- function(exp, kind="cluster") {
  if (kind=='subcluster') {
    fn <- '_pairwiseMarkers.RDS'
  } else {
    fn <- "_FirstRound_pairwiseMarkers.RDS"
  }
  pairwiseMarkers <- readRDS(sprintf("%s/%s%s",exp$exp.dir,exp$base,fn))
  
  as_tibble(mutate(ldply(pairwiseMarkers, function(cluster) {
    mutate(ldply(cluster, function(df) {
      df
    }, .id='other.cluster'), other.cluster=as.character(other.cluster))
  }, .id="cluster"), cluster=as.character(cluster), exp.label=exp$exp.label))
}

lapply(c('cluster','subcluster'), function(kind) {
  dlply(experiments, .(exp.label), function(exp) {
    log(glue("Reading {exp$exp.label} {kind} pairwise Markers"))
    markers <- read.pairwiseMarkers(exp,kind)
    all.genes <- unique(c(all.genes, markers$GENE))
    
    if (kind=='cluster') {
      markers <- mutate(markers, cluster=factor(cluster,levels=levels(cluster.names_$cluster)),
                        other.cluster=factor(other.cluster,levels=levels(cluster.names_$cluster)))
    } else {
      markers <- mutate(markers, cluster=factor(cluster,levels=levels(subcluster.names_$subcluster)),
                        other.cluster=factor(other.cluster,levels=levels(subcluster.names_$subcluster)))
    }
    
    # save each comparison to a separate file
    ddply(markers, .(cluster), function(df) {
      fn.dir <- glue("{markers.dir}/{exp$exp.label}")
      suppressWarnings(dir.create(fn.dir))
      fn <- glue("{fn.dir}/{first(df$cluster)}.{kind}.pairwise.markers.RDS")
      log(glue("Writing {fn}"))
      saveRDS(df, fn)
    })
  }) 
  NULL
})

saveRDS(all.genes, glue("{markers.dir}/all.genes.RDS"))


saveRDS(filter(subclusterMarkers, pval < 1e-200)$GENE %>% union(filter(clusterMarkers, pval < 1e-200)$GENE), glue("{markers.dir}/pval-200.genes.RDS"))

