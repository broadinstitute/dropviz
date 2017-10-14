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
    
    if (!is.null(progress)) progress$set(value=0.6, detail=glue("Fold ratios, p-vals and conf ints for {nrow(x)} genes"))
    
    target.total <- sum(x$target.sum)
    
    ## binom tests are subtly different. The pval is calculated based
    ## on how far target.sum is from
    ## sum(target.sum)/(sum(target.sum)+sum(comparison.sum)), i.e. the
    ## expected random counts based on the size of the two sets if
    ## there was no difference in expression.
    ##
    ## But the L and R confidence is based on the range of probable
    ## counts for target.sum IF target.sum/(target.sum+comparison.sum)
    ## IS the true proportion probability.
    ##
    ## E.g., there may be 10000 transcripts in target (A) and 20000 in
    ## target (B). For gene G, we might observe 10 transcripts in A
    ## and 30 in B for a total of 40. Based on the background counts
    ## (10000 vs 20000), we expect 10000/30000=1/3 of counts in A and
    ## 2/3 in B. So for gene G, we'd expect 1/3*40=13. 10 is close to
    ## 13, so binom.test(10, 40, 1/3)$p.value => 0.316
    ##
    ## But suppose that G is differentially expressed in A vs B. Then
    ## the proportion is estimated at 10/40 = .25. If that's the case
    ## then there is a 95% confidence we would observe counts between
    ## qbinom(0.025, 40, 1/4) and qbinom(0.975, 40, 1/4), i.e. [5,16].
    ## This latter confidence interval measure does not seem
    ## interesting to me because it's merely a function of the total
    ## transcript count, so confidence intervals are tighter as the
    ## observed counts increase. 
    ## 

    x <- mutate(x, 
                log.target.u=log.transform(target.u), 
                log.comparison.u=log.transform(comparison.u),
                pval=edgeR::binomTest(target.sum, comparison.sum),
                fc=log.target.u-log.comparison.u, 
                fc.disp=exp(fc),
                target.sum.L=qbinom(0.025, target.sum+comparison.sum, target.sum/(target.sum+comparison.sum)),
                target.sum.R=qbinom(0.975, target.sum+comparison.sum, target.sum/(target.sum++comparison.sum)),
                target.sum.per.10k=per.10k(target.sum/target.total),
                target.sum.L.per.10k=per.10k(target.sum.L/target.total),
                target.sum.R.per.10k=per.10k(target.sum.R/target.total))

    if (!is.null(progress)) progress$set(0.8, detail=glue("Cacheing pairwise data"))
    
    saveRDS(x, file=cache.file)
  }
  x
}
