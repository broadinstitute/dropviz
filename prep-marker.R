# combine all the pairwise cx-Ncx pair data 
# must run prep-metacells.R before running this

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

all.pairs <- function(exp, cx.fn, lvls) {
  ldply(cx.fn, function(fn) {
    # read cx from fn
    m <- glue("^{exp$exp.abbrev}\\.([^\\.]+)")
    cx <- factor(regmatches(fn, regexec(m, fn))[[1]][2], levels=lvls)
    
    readRDS(glue("{pairs.dir}/{fn}")) %>% mutate(cx=cx, exp.label=exp$exp.label)
    
  }) %>% as_tibble
}

# store all of the comparison against region in single files and save all gene names
all.genes <- 
dlply(experiments, .(exp.label), function(exp) {
  pairs.fn <- list.files(pairs.dir,glue("^{exp$exp.abbrev}.*\\.N.*\\.RDS"))
  subclusters <- grep('-', pairs.fn, value=TRUE)  # hack: subclusters have dashes in the filename
  clusters <- grep('-', pairs.fn, value=TRUE, invert = TRUE)
  
  all.cluster.pairs <- all.pairs(exp, clusters, levels(cluster.names_$cluster))
  write.log(glue("{nrow(all.cluster.pairs)} pairs, {length(unique(all.cluster.pairs$gene))} genes"))
  saveRDS(all.cluster.pairs, glue("{meta.dir}/{exp$exp.label}.gene.clusters.RDS"))
  all.subcluster.pairs <- all.pairs(exp, subclusters, levels(subcluster.names_$subcluster))
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
      base::saveRDS(., fn)
      tibble(fn=as.character(fn))
    })
    
  })
})

saveRDS(uniq.genes, glue("{markers.dir}/expressed.genes.RDS"))

# retrieve comprehensive gene symbols used
all.genes <- ddply(experiments, .(exp.label), function(exp) {
  raw.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.raw.dge.RDS"))
  raw <- readRDS(raw.fn) # rows genes, cols cells in this major cluster
  df <- data.frame(gene=rownames(raw))
  rm(raw)
  df
}) %>% unique
saveRDS(all.genes$gene, glue("{markers.dir}/all.genes.RDS"))

top.genes <-
  ldply(c('clusters','subclusters'), function(kind) {
    ddply(experiments, .(exp.label), function(exp) {
      filter(readRDS(glue("{meta.dir}/{exp$exp.label}.gene.{kind}.RDS")),
             pval<1e-200 & fc.disp > 2 & !grepl('Rik\\d?$',gene,perl = TRUE)) %>% select(gene)
    })
  })

top.genes.symbols <- unique(top.genes$gene)
write.log(glue("Writing {length(top.genes.symbols)} most significantly differentially expressed to {markers.dir}/top_genes.RDS"))
saveRDS(top.genes.symbols, file=glue("{markers.dir}/top_genes.RDS"))

## environment - tolower(symbol) => [ Symbol, Description ]

## retrieved from /broad/mccarroll/software/metadata/individual_reference/mm10/mm10.gene_descriptions.txt via James Nemesh <nemesh@broadinstitute.org>
gene.descriptions <- read.delim("data/gene_descriptions.txt.gz") %>% mutate(Description=sub(' \\[.*','',Description))
write.log(glue("Read {nrow(gene.descriptions)} gene descriptions"))
gene.descriptions <- inner_join(gene.descriptions, all.genes, by=c(Associated.Gene.Name='gene')) # limit to genes used in analysis
write.log(glue("Reduced to {nrow(gene.descriptions)} gene descriptions"))
gene.list <- setNames(lapply(1:nrow(gene.descriptions), function(i) c(gene.descriptions$Associated.Gene.Name[i], gene.descriptions$Description[i])), tolower(gene.descriptions$Associated.Gene.Name))
gene.dict <- list2env(gene.list)
saveRDS(gene.dict, file=glue("{prep.dir}/markers/gene.dict.RDS"))

