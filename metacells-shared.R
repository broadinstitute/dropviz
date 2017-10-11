# metacell functions that are used both online and offline

suppressWarnings(dir.create(glue("{cache.dir}/metacells"), recursive = TRUE))

per.10k <- function(x) 10000*x+1
log.transform <- function(x) log(per.10k(x))

# User can specify >1 cx or cmp.cx. generate weighted means and sum sums.
merge.cxs <- function(cxs, kind) {

  # data is stored per exp.label. Read once for all clusters in that region.
  means.sums <-
    dlply(cxs, .(exp.label), function(cxs.exp) {
      fn.means <- glue("{prep.dir}/metacells/{first(cxs.exp$exp.label)}.{kind}.means.RDS")
      fn.sums <- glue("{prep.dir}/metacells/{first(cxs.exp$exp.label)}.{kind}.sums.RDS")
      
      means <- select(readRDS(fn.means), c('gene',cxs.exp$cx))
      
      sums <- select(readRDS(fn.sums), c('gene', cxs.exp$cx))
      
      list(means=means, sums=sums)
    })
  
  # combine into a single table with only common genes
  inner_join_by_gene <- function(x,y) inner_join(x,y,by='gene')
  means <- Reduce(inner_join_by_gene, lapply(means.sums, function(ms) ms$means))
  sums <- Reduce(inner_join_by_gene, lapply(means.sums, function(ms) ms$sums))
    
  # get relative sizes of clusters
  totals <- apply(sums[2:ncol(sums)], 2, sum)
  grand.total <- sum(totals)
    
  means.vals <- apply(means[2:ncol(means)], 1, function(row) sum(row*totals)/grand.total)
  sums.vals <- apply(sums[2:ncol(sums)], 1, sum)
    
  list(means=tibble(gene=means$gene, means=means.vals), sums=tibble(gene=sums$gene, sums=sums.vals))
}

compute.pair <- function(exp.label, cx, cmp.exp.label, cmp.cx, kind, use.cached=TRUE) {
  
  targets <- tibble(exp.label=exp.label, cx=cx)
  comparisons <- tibble(exp.label=cmp.exp.label, cx=cmp.cx)

  target.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label%in%targets$exp.label]}.{targets$cx}"),collapse='+')
  comparison.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label%in%comparisons$exp.label]}.{comparisons$cx}"),collapse='+')
  
  cache.file <- glue("{cache.dir}/metacells/{target.names}.vs.{comparison.names}.RDS")
  
  if (use.cached && file.exists(cache.file)) {
    x <- readRDS(cache.file)
  } else {
    progress <- shiny.progress(glue("{kind} pairwise - {target.names} vs {comparison.names}"))
    if (!is.null(progress)) on.exit(progress$close())
    
    write.log(glue("Computing pairwise {target.names} vs {comparison.names}"))
    
    if (!is.null(progress)) {
      progress$inc(0.3, detail=glue("Reading means and sums from disk"))
    }

    means.sums.tgt <- merge.cxs(targets, kind)
    means.tgt <- means.sums.tgt$means
    sums.tgt <- means.sums.tgt$sums
    
    means.sums.cmp <- merge.cxs(comparisons, kind)
    means.cmp <- means.sums.cmp$means
    sums.cmp <- means.sums.cmp$sums

    common.genes <- intersect(means.tgt$gene, means.cmp$gene)
    means.tgt <- filter(means.tgt, gene %in% common.genes)
    means.cmp <- filter(means.cmp, gene %in% common.genes)
    sums.tgt <- filter(sums.tgt, gene %in% common.genes)
    sums.cmp <- filter(sums.cmp, gene %in% common.genes)
      
    x <- inner_join(
      inner_join(means.tgt, means.cmp, by='gene') %>% setNames(c('gene','target.u','comparison.u')),
      inner_join(sums.tgt, sums.cmp, by='gene') %>% setNames(c('gene','target.sum','comparison.sum')), by='gene')
    
    if (!is.null(progress)) progress$set(value=0.6, detail=glue("Fold ratio and p-vals for {nrow(x)} genes"))
    
    x <- mutate(x, 
                log.target.u=log.transform(target.u), 
                log.comparison.u=log.transform(comparison.u),
                pval=edgeR::binomTest(target.sum, comparison.sum),
                fc=log.target.u-log.comparison.u, 
                fc.disp=exp(fc),
                size=round(per.10k(target.u)+per.10k(comparison.u)),
                prob=per.10k(target.u)/size,
                log.target.u.L=log(qbinom(0.025, size, prob)),
                log.target.u.R=log(qbinom(0.975, size, prob)))

    if (!is.null(progress)) progress$set(0.8, detail=glue("Cacheing pairwise data"))
    
    saveRDS(x, file=cache.file)
  }
  x
}
