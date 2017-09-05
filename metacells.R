suppressWarnings(dir.create(glue("{cache.dir}/metacells"), recursive = TRUE))


means.and.sums <- function(fn.means, fn.sums, cx, cmp.cx, cache.file, progress) {
  if (file.exists(cache.file)) {
    x <- readRDS(cache.file)
  } else {
    write.log(glue("Computing pairwise {cx} vs {cmp.cx}"))
    progress$set(value=0.3, message=glue("{cx} vs {cmp.cx}"), detail=glue("Reading means and sums from disk"))
    x <- bind_cols(readRDS(fn.means)[,c('gene',cx,cmp.cx)] %>% setNames(c('gene','target.u','comparison.u')),
                   readRDS(fn.sums)[,c(cx,cmp.cx)] %>% setNames(c('target.sum','comparison.sum')))
    progress$set(value=0.6, detail=glue("Fold ratio and p-vals for {nrow(x)} genes"))
    x <- mutate(x, log.target.u=log(10000*target.u+1), 
                log.comparison.u=log(10000*comparison.u+1),
                pval=edgeR::binomTest(target.sum, comparison.sum),
                fc=log.target.u-log.comparison.u, 
                fc.disp=exp(fc))
    progress$set(value=0.8, detail=glue("Cacheing pairwise data"))
    saveRDS(x, file=cache.file)
  }

  mutate(x, fc.thresh=( if (input$expr.filter.opt %in% 'fc') { (abs(fc) > log(input$fold.change)) } else { fc > log(input$fold.change) } ),
         pval.thresh=pval < 10^input$pval.thresh,
         amt.thresh=( (log.target.u > input$min.amt.within) & (log.comparison.u < input$max.amt.without) ),
         user.selected=gene %in% input$user.genes,
         expr.pass = (if (input$expr.filter.opt == 'either') { (fc.thresh & pval.thresh) | amt.thresh } else {
           if (input$expr.filter.opt=='both') { fc.thresh & pval.thresh & amt.thresh } else {
             if (input$expr.filter.opt == 'fc') { fc.thresh & pval.thresh } else { amt.thresh } 
           }
         } )
  )
}

# return a tibble of gene expression (region.disp, gene, target, comparison) 
# for all the cells in the current.cluster() and comparison.cluster() 
# if the comparison is global, then expression is all cells except the target in the target region.
cluster.metacells.selected <- reactive({
log.reactive("fn: cluster.metacells.selected")
  req(input$current.cluster)
  
  fn.means <- glue("{prep.dir}/metacells/{current.cluster()$exp.label}.cluster.means.RDS")
  fn.sums <- glue("{prep.dir}/metacells/{current.cluster()$exp.label}.cluster.sums.RDS")

  cx <- as.character(current.cluster()$cluster)
  cmp.cx <- ifelse(comparison.cluster()$cluster == 'global', paste0('N', current.cluster()$cluster), as.character(comparison.cluster()$cluster))
  cache.file <- glue("{cache.dir}/metacells/{current.cluster()$exp.label}.{cx}.{cmp.cx}.RDS")
  region.disp <- current.cluster()$region.disp

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  means.and.sums(fn.means, fn.sums, cx, cmp.cx, cache.file, progress) %>% mutate(region.disp=region.disp)  
})

# same as above, but for subclusters
subcluster.metacells.selected <- reactive({
log.reactive("fn: subcluster.metacells.selected")
  req(input$current.subcluster)

  fn.means <- glue("{prep.dir}/metacells/{current.subcluster()$exp.label}.subcluster.means.RDS")
  fn.sums <- glue("{prep.dir}/metacells/{current.subcluster()$exp.label}.subcluster.sums.RDS")
  
  cx <- as.character(current.subcluster()$subcluster)
  cmp.cx <- ifelse(comparison.subcluster()$subcluster == 'global', paste0('N', current.subcluster()$subcluster), as.character(comparison.subcluster()$subcluster))
  cache.file <- glue("{cache.dir}/metacells/{current.subcluster()$exp.label}.{cx}.{cmp.cx}.RDS")
  region.disp <- current.subcluster()$region.disp
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  means.and.sums(fn.means, fn.sums, cx, cmp.cx, cache.file, progress) %>% mutate(region.disp=region.disp)
})

# returns a closure of a scatter plot of target vs comparison
scatter.plot <- function(target, comparison, cx.metacells, cx.markers, cx.selected.markers, kind, return.closure=FALSE) {
  function() {
    target.cx <- target[[kind]]
    target.label <- filter(cx.names(), exp.label==target$exp.label, cx==as.character(target.cx))$cx.disp
    
    compare.cx <- comparison[[kind]]
    if (compare.cx == 'global') {
      compare.cx <- paste0('N',target.cx)
      compare.label <- filter(region.names(), exp.label==target$exp.label)$region.disp
    } else {
      compare.label <- filter(cx.names(), exp.label==comparison$exp.label, cx==as.character(compare.cx))$cx.disp
    }
    
    cell.pairs <- mutate(cx.metacells, 
                         pval=ifelse(pval==0, .Machine$double.xmin, pval))  # to allow log

    cell.pairs$pval.disp <- log(cell.pairs$pval)
    cell.pairs$criteria.thresh <- cell.pairs$gene %in% cx.markers$gene
    cell.pairs$gene.tbl.selected <- cell.pairs$gene %in% cx.selected.markers$gene
    cell.pairs$kind <- factor(with(cell.pairs, ifelse(user.selected, 'User Set', ifelse(gene.tbl.selected, 'Selected', ifelse(criteria.thresh, 'Pass Criteria', 'Other')))), levels=c('Other','Pass Criteria','Selected','User Set'))
    cell.pairs$size <- factor(ifelse(cell.pairs$kind %in% c('Selected','User Set','Pass Criteria'),'Filtered','Other'), levels=c('Other','Filtered'))
    cell.pairs$gene.label <- ifelse(cell.pairs$gene.tbl.selected | cell.pairs$user.selected, cell.pairs$gene, NA_character_)
    
    opt.scatter.gene.labels <- input$opt.scatter.gene.labels
    opt.expr.size <- input$opt.expr.size
    p.func <- function() {
      require(ggplot2); require(dplyr)
      
      gene.labels <- (
        if (opt.scatter.gene.labels && any(!is.na(cell.pairs$gene.label))) {
          geom_label(aes(x=log.target.u+.3, label=gene.label), size=3, show.legend = FALSE)
        } else {
          geom_blank()
        }
      )
      
      if (input$expr.filter.opt %in% c('both','fc','either')) {
        if (input$expr.filter.opt == 'fc') {
          fc.line1 <- geom_abline(intercept=log(input$fold.change), slope=1) 
        } else {
          fc.line1 <- geom_blank()
        }
        fc.line2 <- geom_abline(intercept=-log(input$fold.change), slope=1)
      } else {
        fc.line1 <- fc.line2 <- geom_blank()
      }
      
      if (input$expr.filter.opt %in% c('both','amt','either')) {
        amt.line1 <- geom_hline(yintercept = input$max.amt.without) 
        amt.line2 <- geom_vline(xintercept = input$min.amt.within)
      } else {
        amt.line1 <- amt.line2 <- geom_blank()
      }

      plot.gg <- ggplot(cell.pairs, aes(x=log.target.u, y=log.comparison.u, color=pval.disp, size=size)) + geom_point() + 
        fc.line1 + fc.line2 +
        amt.line1 + amt.line2 + 
        gene.labels +
        geom_point(data=filter(cell.pairs, !is.na(gene.label) & expr.pass), color='green') +
        geom_point(data=filter(cell.pairs, !is.na(gene.label) & !expr.pass), color='red') +
        scale_color_continuous(name='p-val exp') + 
        scale_size_manual(guide="none",values=c('Other'=.5,'Filtered'=opt.expr.size)) + 
        facet_wrap(~region.disp) + xlab(paste(target.label,"\n(log normal mean)")) + ylab(paste(compare.label,"\n(log normal mean)"))

    }
    
    if (return.closure) {
      p.func
    } else {
      p.func()
    }
    
  }
    
}

# draw a scatter plot of expression levels between the target and comparison set
# All genes from the differentially expressed analysis are displayed.
# In addition, transcript counts for any user genes are retrieved from file.
output$gene.expr.scatter.cluster <- renderImage({
  if (isTruthy(input$current.cluster) && isTruthy(input$comparison.cluster)) {
    cluster.scatter.plot <- scatter.plot(current.cluster(), comparison.cluster(), cluster.metacells.selected(), cluster.markers(), cluster.markers.selected(), 'cluster')
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.cluster()$exp.label, 
                        current.cluster()$cluster, comparison.cluster()$cluster, cluster.markers()$gene, 
                        cluster.markers.selected()$gene, input$expr.filter.opt, input$fold.change,
                        input$max.amt.without, input$min.amt.within))
  } else {
    cluster.scatter.plot <- function() plot.text("Choose a Target Cluster and Comparison")    
    key.str <- 'missing_clusters'
  }
  
  height <- img.size.round(session$clientData[[glue("output_gene.expr.scatter.cluster_height")]])
  width <- height
  
  renderCacheImage(cluster.scatter.plot, glue("scatter_cluster_{key.str}"), width, height)
}, deleteFile=FALSE)

output$gene.expr.scatter.cluster.dl <- downloadHandler(filename="scatter.zip", 
                                                       content= function(file) {
                                                         req(isTruthy(input$current.cluster) && isTruthy(input$comparison.cluster))
                                                         cluster.scatter.plot <- scatter.plot(current.cluster(), comparison.cluster(), cluster.metacells.selected(), cluster.markers(), cluster.markers.selected(), 'cluster', return.closure = TRUE)()
                                                         send.zip(cluster.scatter.plot, 'scatter', file)
                                                       })

output$gene.expr.scatter.subcluster <- renderImage({
  if (isTruthy(input$current.subcluster) && isTruthy(input$comparison.subcluster)) {
    subcluster.scatter.plot <- scatter.plot(current.subcluster(), comparison.subcluster(), subcluster.metacells.selected(), subcluster.markers(), subcluster.markers.selected(), 'subcluster')
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.subcluster()$exp.label, current.subcluster()$subcluster, comparison.subcluster()$subcluster, subcluster.markers()$gene, subcluster.markers.selected()$gene))
  } else {
    subcluster.scatter.plot <- function() plot.text("Choose a Target Subcluster and Comparison")    
    key.str <- 'missing_subclusters'
  }
  
  height <- img.size.round(session$clientData[[glue("output_gene.expr.scatter.subcluster_height")]])
  width <- height
  
  renderCacheImage(subcluster.scatter.plot, glue("scatter_subcluster_{key.str}"), width, height)
}, deleteFile=FALSE)

output$gene.expr.scatter.subcluster.dl <- downloadHandler(filename="scatter.zip", 
                                                          content= function(file) {
                                                            req(isTruthy(input$current.subcluster) && isTruthy(input$comparison.subcluster))
                                                            subcluster.scatter.plot <- scatter.plot(current.subcluster(), comparison.subcluster(), subcluster.metacells.selected(), subcluster.markers(), subcluster.markers.selected(), 'subcluster', return.closure = TRUE)()
                                                            send.zip(subcluster.scatter.plot, 'scatter', file)
                                                          })
