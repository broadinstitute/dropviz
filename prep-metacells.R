# create new metacells files that include transcript counts for all clusters except one, labeled "N#".

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

ddply(experiments, .(exp.label), function(exp) {
  log("Generating cluster metacell data for  ",exp$exp.label)
  mc <- readRDS(glue("{exp$exp.dir}/metacells/{exp$base}.metacells.RDS"))
  mc.sum <- apply(mc, 1, sum)
  mc <- bind_cols(tibble(gene=rownames(mc)),as_tibble(mc))
  mc <- mk.neg.cols(mc, mc.sum)
  saveRDS(mc, file=glue("{meta.dir}/{exp$exp.label}.cluster.metacells.RDS"))
  
  log("Generating subcluster metacell data for  ",exp$exp.label)
  mc <- readRDS(glue("{exp$exp.dir}/metacells/{exp$base}.subcluster.metacells.RDS"))
  mc.sum <- apply(mc, 1, sum)
  mc <- bind_cols(tibble(gene=rownames(mc)),as_tibble(mc))
  mc <- mk.neg.cols(mc, mc.sum)
  saveRDS(mc, file=glue("{meta.dir}/{exp$exp.label}.subcluster.metacells.RDS"))
})
