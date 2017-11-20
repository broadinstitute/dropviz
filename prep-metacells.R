# Create new "metacell" structures: 
# - the sum of transcript counts for all genes in group,  "metacell.sums"
# - the mean of normalized counts for all genes in group (where counts are normalized to sum to 1 within cell)  metacell.means
# where group is each cluster and !cluster

# Also pre-compute all pairwise stats for each cluster/subcluster vs
# region. This is used by prep-marker to generate fast lookup of
# differentially expressed genes and also provides quick app response
# for diffex table because default comparison is the region. The code
# that is used is the same code that's used during runtime for
# arbitrary cx1 vs cx2 computations. These cx vs Ncx calcs are stored
# in prep.dir by locally reassigning cache.dir here and in
# prep-marker.R. Later those pre-computed values will be transferred to
# local cache.dir for each new server instance.

# The calculations here are (hopefully) consistent with the statistics used in IcaCluster_Functions.R:diffexp

# scaled DGE: total expression f each cell sum to 1.  Genes X Cells
# raw DGE: Genes X Cells (raw transcript counts)
# metacell: sum of raw transcript counts per cluster - Genes x Clusters - same as raw.row.sums in compute_diffexp_mean_and_sum
# metacell.means: mean of scaled (normalized) transcripts per cluster - Genes x Cluster - same as norm.row.means in compute_diffexp_mean_and_sum

require(Matrix)
source("global.R")
source("metacells-shared.R")

meta.dir <- glue('{prep.dir}/metacells')
suppressWarnings(dir.create(meta.dir, recursive = TRUE))

pairs.dir <- glue("{prep.dir}/pairs")
suppressWarnings(dir.create(pairs.dir, recursive = TRUE))

mk.neg.cols <- function(mc, mc.sum) {
  cols <- 2:ncol(mc)
  lapply(cols, function(i) {
    new.col.name <- paste0('N',colnames(mc)[i])
    mc[[new.col.name]] <<- mc.sum - mc[[i]]
  })
  mc
}

# adds "N{group}" labels as all the cells except the group
mk.neg.group <- function(group.assignment) {
  ldply(levels(group.assignment$group), function(group) {
    data.frame(cell=group.assignment$cell[group.assignment$group != group], group=paste0('N',as.character(group)))
  })
}

# returns table of rows are gene sums and columns are group (cluster, subcluster, !cluster, !subcluster)
compute.dge.stat <- function(group.assignment, mtx, func) {
  
  fns <- 
    dlply(group.assignment, .(group), function(df, func) {
      write.log(first(df$group))
      fn <- tempfile()
      x <- func(mtx[,df$cell, drop=FALSE])
      saveRDS(x, file=fn)
      rm(x)
      gc()
      fn
    }, func=func)
  
  bind_cols(tibble(gene=rownames(mtx)),
            lapply(fns, function(fn) {
              readRDS(fn)
            }) %>% as_tibble
  )
  
}



ddply(experiments, .(exp.label), function(exp) {
  flag <- glue("{meta.dir}/{exp$exp.label}.flag")
  if (!file.exists(flag)) {
    close(file(flag,"w"))
    write.log("Generating cluster metacell data for ",exp$exp.label)
    # convert factor to table
    cluster.cell.assign <- readRDS(sprintf("%s/assign/%s.cluster.assign.RDS", exp$exp.dir, exp$base))
    cluster.cell.assign.tbl <- na.omit(tibble(cell=names(cluster.cell.assign), group=factor(cluster.cell.assign, levels=levels(cell.types$cluster))))
    
    saveRDS(group_by(cluster.cell.assign.tbl, group) %>% summarize(count=length(cell)), sprintf(glue("{meta.dir}/{exp$exp.label}.cluster.counts.RDS")))
    
    cluster.cell.assign.tbl <- rbind(cluster.cell.assign.tbl,
                                     mk.neg.group(cluster.cell.assign.tbl))
    
    subcluster.cell.assign <- readRDS(sprintf("%s/assign/%s.subcluster.assign.RDS", exp$exp.dir, exp$base))
    subcluster.cell.assign.tbl <- na.omit(tibble(cell=names(subcluster.cell.assign), group=factor(subcluster.cell.assign, levels=levels(cell.types$subcluster))))
    
    saveRDS(group_by(subcluster.cell.assign.tbl, group) %>% summarize(count=length(cell)), sprintf(glue("{meta.dir}/{exp$exp.label}.subcluster.counts.RDS")))
    
    subcluster.cell.assign.tbl <- rbind(subcluster.cell.assign.tbl,
                                        mk.neg.group(subcluster.cell.assign.tbl))
    
    scaled.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.scaled.dge.RDS"))
    
    cluster.means.fn <- glue("{meta.dir}/{exp$exp.label}.cluster.means.RDS")
    if (!file.exists(cluster.means.fn) || !getOption("dropviz.prep.cache", default = TRUE)) {
      scaled <- readRDS(scaled.fn) # rows genes, cols cells in this major cluster
      cluster.means <- compute.dge.stat(cluster.cell.assign.tbl, scaled, rowMeans)
      saveRDS(cluster.means, file=cluster.means.fn)
      rm(scaled, cluster.means)
    }
    subcluster.means.fn <- glue("{meta.dir}/{exp$exp.label}.subcluster.means.RDS")
    if (!file.exists(subcluster.means.fn) || !getOption("dropviz.prep.cache", default = TRUE)) {
      scaled <- readRDS(scaled.fn) # rows genes, cols cells in this major cluster
      subcluster.means <- compute.dge.stat(subcluster.cell.assign.tbl, scaled, rowMeans)
      saveRDS(subcluster.means, file=subcluster.means.fn)
      rm(scaled, subcluster.means)
    }

    
    raw.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.raw.dge.RDS"))
    
    cluster.sums.fn <- glue("{meta.dir}/{exp$exp.label}.cluster.sums.RDS")
    if (!file.exists(cluster.sums.fn) || !getOption("dropviz.prep.cache", default = TRUE)) {
      raw <- readRDS(raw.fn) # rows genes, cols cells in this major cluster
      cluster.sums <- compute.dge.stat(cluster.cell.assign.tbl, raw, rowSums)
      saveRDS(cluster.sums, file=cluster.sums.fn)
      rm(raw, cluster.sums)
    }
    subcluster.sums.fn <- glue("{meta.dir}/{exp$exp.label}.subcluster.sums.RDS")
    if (!file.exists(subcluster.sums.fn) || !getOption("dropviz.prep.cache", default = TRUE)) {
      raw <- readRDS(raw.fn) # rows genes, cols cells in this major cluster
      subcluster.sums <- compute.dge.stat(subcluster.cell.assign.tbl, raw, rowSums)
      saveRDS(subcluster.sums, file=subcluster.sums.fn)
      rm(raw, subcluster.sums)
    }
    
    gc()
    
    unlink(flag)
  }
})

do.pairwise.Ncx <- function(exp.label, groups, kind) {
  lapply(groups, function(cx) {
    cmp.cx <- paste0('N',cx)
    compute.pair(exp.label, as.character(cx), exp.label, cmp.cx, kind, use.cached = getOption("dropviz.prep.cache", default = TRUE), pairs.dir=pairs.dir)
  })
}

ddply(experiments, .(exp.label), function(exp) {
  do.pairwise.Ncx(exp$exp.label, unique(cell.types$cluster[cell.types$exp.label==exp$exp.label]),'cluster')
  do.pairwise.Ncx(exp$exp.label, cell.types$subcluster[cell.types$exp.label==exp$exp.label],'subcluster')
  exp
})
