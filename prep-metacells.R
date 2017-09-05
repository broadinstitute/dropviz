# Create new "metacell" structures: 
# - the sum of transcript counts for all genes in group,  "metacell.sums"
# - the mean of normalized counts for all genes in group (where counts are normalized to sum to 1 within cell)  metacell.means
# where group is each cluster and !cluster

# This is (hopefully) consistent with the statistics used in IcaCluster_Functions.R:diffexp

# scaled DGE: total expression f each cell sum to 1.  Genes X Cells
# raw DGE: Genes X Cells (raw transcript counts)
# metacell: sum of raw transcript counts per cluster - Genes x Clusters - same as raw.row.sums in compute_diffexp_mean_and_sum
# metacell.means: mean of scaled (normalized) transcripts per cluster - Genes x Cluster - same as norm.row.means in compute_diffexp_mean_and_sum

source("global.R")

meta.dir <- glue('{prep.dir}/metacells')
suppressWarnings(dir.create(meta.dir, recursive = TRUE))

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
  
  # due to RAM limits, write each to disk and then combine afterwards
  fns <- 
    dlply(group.assignment, .(group), function(df, func) {
      write.log(as.character(substitute(func))," group ",first(df$group))
      fn <- tempfile()
      saveRDS(func(mtx[,df$cell, drop=FALSE]), file=fn)
      fn
    }, func=func)
  
  bind_cols(tibble(gene=rownames(mtx)),
            lapply(fns, function(fn) {
              readRDS(fn)
            }) %>% as_tibble
  )
  
}



ddply(experiments, .(exp.label), function(exp) {
  write.log("Generating cluster metacell data for  ",exp$exp.label)
  # convert factor to table
  cluster.cell.assign <- readRDS(sprintf("%s/assign/%s.cluster.assign.RDS", exp$exp.dir, exp$base))
  cluster.cell.assign.tbl <- na.omit(tibble(cell=names(cluster.cell.assign), group=factor(cluster.cell.assign, levels=levels(cell.types$cluster))))
  cluster.cell.assign.tbl <- rbind(cluster.cell.assign.tbl,
                                   mk.neg.group(cluster.cell.assign.tbl))
  
  subcluster.cell.assign <- readRDS(sprintf("%s/assign/%s.subcluster.assign.RDS", exp$exp.dir, exp$base))
  subcluster.cell.assign.tbl <- na.omit(tibble(cell=names(subcluster.cell.assign), group=factor(subcluster.cell.assign, levels=levels(cell.types$subcluster))))
  subcluster.cell.assign.tbl <- rbind(subcluster.cell.assign.tbl,
                                      mk.neg.group(subcluster.cell.assign.tbl))
  
  scaled.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.scaled.dge.RDS"))
  scaled <- readRDS(scaled.fn) # rows genes, cols cells in this major cluster

  cluster.means <- compute.dge.stat(cluster.cell.assign.tbl, scaled, rowMeans)
  saveRDS(cluster.means, file=glue("{meta.dir}/{exp$exp.label}.cluster.means.RDS"))
  subcluster.means <- compute.dge.stat(subcluster.cell.assign.tbl, scaled, rowMeans)
  saveRDS(subcluster.means, file=glue("{meta.dir}/{exp$exp.label}.subcluster.means.RDS"))
  rm(scaled)
  
  raw.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.raw.dge.RDS"))
  raw <- readRDS(raw.fn) # rows genes, cols cells in this major cluster
  
  cluster.sums <- compute.dge.stat(cluster.cell.assign.tbl, raw, rowSums)
  saveRDS(cluster.sums, file=glue("{meta.dir}/{exp$exp.label}.cluster.sums.RDS"))
  subcluster.sums <- compute.dge.stat(subcluster.cell.assign.tbl, raw, rowSums)
  saveRDS(subcluster.sums, file=glue("{meta.dir}/{exp$exp.label}.subcluster.sums.RDS"))
  rm(raw)
  
  gc()
  
})
