# metacell functions that are used both online and offline

suppressWarnings(dir.create(glue("{cache.dir}/metacells"), recursive = TRUE))



compute.pair <- function(exp.label, cx, cmp.cx, kind, progress=NULL) {
  fn.means <- glue("{prep.dir}/metacells/{exp.label}.{kind}.means.RDS")
  fn.sums <- glue("{prep.dir}/metacells/{exp.label}.{kind}.sums.RDS")
  cache.file <- glue("{cache.dir}/metacells/{exp.label}.{cx}.{cmp.cx}.RDS")
  
  if (file.exists(cache.file)) {
    x <- readRDS(cache.file)
  } else {
    write.log(glue("Computing pairwise {cx} vs {cmp.cx}"))
    
    if (!is.null(progress)) {
      write.log("Set progress")
      progress$set(value=0.3, message=glue("{cx} vs {cmp.cx}"), detail=glue("Reading means and sums from disk"))
    }

    means <- readRDS(fn.means)
    sums <- readRDS(fn.sums)
    if (is.null(means[[cx]]) || is.null(means[[cmp.cx]])) {
      write.log("Missing data. Skipping.")
      x <- tibble(gene=character(0), target.u=double(0), comparison.u=double(0), 
                  target.sum=double(0), comparison.sum=double(0),
                  log.target.u=double(0), log.comparison.u=double(0),
                  pval=double(0), fc=double(0), fc.disp=double(0))
    } else {
      x <- bind_cols(means[,c('gene',cx,cmp.cx)] %>% setNames(c('gene','target.u','comparison.u')),
                     sums[,c(cx,cmp.cx)] %>% setNames(c('target.sum','comparison.sum')))
      
      if (!is.null(progress)) progress$set(value=0.6, detail=glue("Fold ratio and p-vals for {nrow(x)} genes"))
      
      x <- mutate(x, log.target.u=log(10000*target.u+1), 
                  log.comparison.u=log(10000*comparison.u+1),
                  pval=edgeR::binomTest(target.sum, comparison.sum),
                  fc=log.target.u-log.comparison.u, 
                  fc.disp=exp(fc))
    }
    
    
    if (!is.null(progress)) progress$set(value=0.8, detail=glue("Cacheing pairwise data"))
    
    saveRDS(x, file=cache.file)
  }
  x
}
