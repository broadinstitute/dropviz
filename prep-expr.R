# save the rows of a DGE file to separate locations for fast access
# <cachedir>/expr/<region>/gene/<gene>.RDS

library(plyr)
library(dplyr)
library(glue)

source("global.R")
all.genes <- readRDS(glue("{prep.dir}/markers/all.genes.RDS"))

expr.dir <- paste0(prep.dir,'/expr')
xy.dir <- paste0(prep.dir,'/tsne')
N <- 500 # chunk size

gene.map <-
  ddply(experiments, .(exp.label), function(exp) {
    write.log(exp$exp.label)
    expr.exp.dir <- glue("{expr.dir}/{exp$exp.label}")
    suppressWarnings(dir.create(expr.exp.dir, recursive = TRUE))
    
    expr.fn <- with(exp, glue("{exp.dir}/dge/{base}.filtered.scaled.dge.RDS"))
    expr <- readRDS(expr.fn) # rows genes, cols cells

    genes.fn <- glue("{expr.exp.dir}/genes.RDS")
    if (!file.exists(genes.fn))
      saveRDS(rownames(expr), genes.fn)
    
    # split the matrix into chunks of N rows at a time, convert to dense, and then write to disk
    x <-
      ldply(seq(1, nrow(expr), by=N), function(i) {
        fn <- glue("{expr.exp.dir}/{i}.RDS")
        if (!file.exists(fn)) {
          write.log(glue("Extracting {N} rows starting at {i}"))
          sub.expr <- expr[i:min(nrow(expr),(i+N-1)),]
          
          write.log(glue("Converting and writing {fn}"))
          saveRDS(as.matrix(sub.expr), fn)
          rm(sub.expr)
          gc()
        }
        
        data.frame(exp.label=exp$exp.label, start.gene=i, end.gene=i+N-1, file=fn)
      })
    rm(expr)
    gc()
    x
  })

saveRDS(gene.map, glue("{expr.dir}/gene.map.RDS"))

# create a tibble of all the non-zero entries
expr.chunks <-
  ddply(gene.map, .(exp.label), function(exp) {
    expr.exp.dir <- glue("{expr.dir}/{first(exp$exp.label)}")
    
    gene.names <- readRDS(glue("{expr.exp.dir}/genes.RDS"))
    
    genes.expr <-  
      ddply(exp, .(start.gene), function(chunk) {
        write.log(glue("Processing {chunk$file}"))  
        out.fn <- glue("{expr.exp.dir}/{chunk$start.gene}_df.RDS")
        
        if (!file.exists(out.fn)) {
          sub.expr <- readRDS(chunk$file)
          chunk.df <-       
            ldply(1:nrow(sub.expr), function(r) {
              non.zero <- sub.expr[r,]!=0
              data.frame(gene=as.factor(rownames(sub.expr)[r]), cell=factor(names(sub.expr[r,non.zero])), transcripts=sub.expr[r,non.zero])
            }) %>% as_tibble
          saveRDS(chunk.df, out.fn)
        }
        data.frame(file=out.fn)
      })
  })

ddply(expr.chunks, .(exp.label), function(exp) {
  out.dir <- glue("{expr.dir}/{first(exp$exp.label)}/gene")
  suppressWarnings(dir.create(out.dir, recursive = TRUE))
  xy <- select(readRDS(glue("{xy.dir}/{first(exp$exp.label)}/global.xy.RDS")), -subcluster)
  
  ddply(exp, .(file), function(chunk) {
    chunk.expr <- readRDS(chunk$file)
    chunk.expr <- filter(chunk.expr, gene %in% all.genes)
    write.log(glue("{chunk$file} has {length(unique(chunk.expr$gene))} filtered genes"))
    ddply(chunk.expr, .(gene), function(g) {
      # replace factors with char. Later we will load and join by character or coerce this to
      # a factor representation using ordinals established elsewhere. Saving the factor with 10,000s of possible values 
      # here wastes a lot of space.
      g <- mutate(g, gene=as.character(gene), cell=as.character(cell)) %>% inner_join(xy, by='cell')
      fn <- glue("{out.dir}/{first(g$gene)}.RDS")
      saveRDS(as_tibble(g), fn)
    })
    
  })
  data.frame()
})

