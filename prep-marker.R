# combine all the pairwise cx-Ncx pair data 

library(plyr)
library(dplyr)
library(glue)

source("global.R")
markers.dir <- glue("{prep.dir}/markers")
genes.dir <- glue("{markers.dir}/genes")
pairs.dir <- glue("{cache.dir}/metacells")
meta.dir <- glue('{prep.dir}/metacells')

suppressWarnings(dir.create(genes.dir, recursive = TRUE))
suppressWarnings(dir.create(glue("{genes.dir}/cluster"), recursive = TRUE))
suppressWarnings(dir.create(glue("{genes.dir}/subcluster"), recursive = TRUE))

all.pairs <- function(exp.label, cx.fn, lvls) {
  ldply(cx.fn, function(fn) {
    # read cx from fn
    m <- glue("^{exp.label}\\.([^\\.]+)")
    cx <- factor(regmatches(fn, regexec(m, fn))[[1]][2], levels=lvls)
    
    readRDS(glue("{pairs.dir}/{fn}")) %>% mutate(cx=cx, exp.label=exp.label)
    
  }) %>% as_tibble
}

# store all of the comparison against region in single files and save all gene names
all.genes <- 
dlply(experiments, .(exp.label), function(exp) {
  pairs.fn <- list.files(pairs.dir,glue("^{exp$exp.label}.*\\.N.*\\.RDS"))
  subclusters <- grep('-', pairs.fn, value=TRUE)  # hack: subclusters have dashes in the filename
  clusters <- grep('-', pairs.fn, value=TRUE, invert = TRUE)
  
  all.cluster.pairs <- all.pairs(exp$exp.label, clusters, levels(cluster.names_$cluster))
  write.log(glue("{nrow(all.cluster.pairs)} pairs, {length(unique(all.cluster.pairs$gene))} genes"))
  saveRDS(all.cluster.pairs, glue("{meta.dir}/{exp$exp.label}.gene.clusters.RDS"))
  all.subcluster.pairs <- all.pairs(exp$exp.label, subclusters, levels(subcluster.names_$subcluster))
  write.log(glue("{nrow(all.subcluster.pairs)} pairs, {length(unique(all.subcluster.pairs$gene))} genes"))
  saveRDS(all.subcluster.pairs, glue("{meta.dir}/{exp$exp.label}.gene.subclusters.RDS"))
  
  unique(all.cluster.pairs$gene)
})
uniq.genes <- unique(unlist(all.genes, use.names=FALSE))

# divide gene names into small groups of 1000 genes
nr <- length(uniq.genes)
n <- 1000
gene.groups <- split(uniq.genes, rep(1:ceiling(nr/n), each=n, length.out=nr))

lapply(gene.groups, function(group) {
  lapply(c('cluster','subcluster'), function(kind) {
    write.log(paste(group[1],group[length(group)],kind))
    gene.set <-
      do.call(bind_rows, 
              lapply(list.files(meta.dir,glue("\\.{kind}s.RDS")), function(fn) filter(readRDS(glue("{meta.dir}/{fn}")), gene %in% group)))
    
    gene.set %>% group_by(gene) %>% do({
      fn <- glue("{genes.dir}/{kind}/{first(.$gene)}.diffexp.RDS")
      saveRDS(., fn)
      tibble(fn=as.character(fn))
    })
    
  })
})

saveRDS(uniq.genes, glue("{markers.dir}/all.genes.RDS"))

#saveRDS(filter(subclusterMarkers, pval < 1e-200)$GENE %>% union(filter(clusterMarkers, pval < 1e-200)$GENE), glue("{markers.dir}/pval-200.genes.RDS"))

