source("metacells-shared.R", local=TRUE)

dir.create(file.path(cache.dir,"metacells"), showWarnings=FALSE)

# adds local attributes to results of compute.pair
cx.pairwise <- function(exp.label, cx, cmp.exp.label, cmp.cx, kind) {

  mutate(compute.pair(exp.label, cx, cmp.exp.label, cmp.cx, kind), 
         fc.thresh=( if (input$expr.filter.opt %in% 'fc') { (abs(fc) > log(input$fold.change)) } else { fc > log(input$fold.change) } ),
         pval.thresh=pval < 10^input$pval.thresh,
         amt.thresh=( (log.target.u > input$min.amt.within) & (log.comparison.u < input$max.amt.without) ),
         user.selected=gene %in% user.manual.genes(),
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
  req(input$current.cluster)
  req(current.cluster()$cluster)
  
  cx <- as.character(current.cluster()$cluster)
  cmp.cx <- ifelse(comparison.cluster()$cluster == 'global', paste0('N', current.cluster()$cluster), as.character(comparison.cluster()$cluster))

  cx.pairwise(current.cluster()$exp.label, cx, comparison.cluster()$exp.label, cmp.cx, 'cluster')
})

# same as above, but for subclusters
subcluster.metacells.selected <- reactive({
  req(input$current.subcluster)
  req(current.subcluster()$subcluster)

  cx <- as.character(current.subcluster()$subcluster)
  cmp.cx <- ifelse(comparison.subcluster()$subcluster == 'global', paste0('N', current.subcluster()$subcluster), as.character(comparison.subcluster()$subcluster))

  cx.pairwise(current.subcluster()$exp.label, cx, comparison.subcluster()$exp.label, cmp.cx, 'subcluster')
})

# returns a closure of a scatter plot of target vs comparison
scatter.plot <- function(target, comparison, cx.metacells, cx.markers, cx.selected.markers, kind, return.closure=FALSE) {
  function(progress=NULL) {
    target$cx <- as.character(target[[kind]])
    target.label <- paste(inner_join(cx.names(), target, by=c('exp.label','cx'))$cx.disp, collapse='+')
    
    comparison$cx <- as.character(comparison[[kind]])
    if (nrow(comparison)==1 && comparison$cx == 'global') {
      comparison$cx <- paste0('N',target$cx)
      compare.label <- filter(region.names(), exp.label==target$exp.label)$region.disp
    } else {
      compare.label <- paste(inner_join(cx.names(), comparison, by=c('exp.label','cx'))$cx.disp, collapse='+')
    }
    
    cell.pairs <- mutate(cx.metacells, 
                         pval=ifelse(pval==0, .Machine$double.xmin, pval),  # to allow taking log
                         target.region=paste(unique(target$region.disp), collapse=' + '),
                         compare.region=paste(unique(comparison$region.disp), collapse=' + ')
    )


    cell.pairs$pval.disp <- log(cell.pairs$pval)
    cell.pairs$criteria.thresh <- cell.pairs$gene %in% cx.markers$gene
    cell.pairs$gene.tbl.selected <- cell.pairs$gene %in% cx.selected.markers$gene
    cell.pairs$kind <- factor(with(cell.pairs, ifelse(user.selected, 'User Set', ifelse(gene.tbl.selected, 'Selected', ifelse(criteria.thresh, 'Pass Criteria', 'Other')))), levels=c('Other','Pass Criteria','Selected','User Set'))
    cell.pairs$size <- factor(ifelse(cell.pairs$kind %in% c('Selected','User Set','Pass Criteria'),'Filtered','Other'), levels=c('Other','Filtered'))
    cell.pairs$gene.label <- ifelse(cell.pairs$gene.tbl.selected | cell.pairs$user.selected, cell.pairs$gene, NA_character_)
    
    opt.scatter.gene.labels <- input$opt.scatter.gene.labels
    opt.expr.size <- input$opt.expr.size
    opt.expr.filter <- input$expr.filter.opt
    opt.fold.change <- input$fold.change
    opt.max.amt.without <- input$max.amt.without
    opt.min.amt.within <- input$min.amt.within
    p.func <- function() {
      source("scatter-plot.R", local=TRUE)
    }
    
    if (return.closure) {
      p.func
    } else {
      if (!is.null(progress)) progress$inc(detail=glue("Rendering {nrow(cell.pairs)} cell pairs in scatter plot"))
      p.func()
    }
    
  }
    
}

# draw a scatter plot of expression levels between the target and comparison set
# All genes from the differentially expressed analysis are displayed.
# In addition, transcript counts for any user genes are retrieved from file.
output$gene.expr.scatter.cluster <- renderImage({
  progress <- shiny.progress('Scatter')
  if (!is.null(progress)) on.exit(progress$close())
  
  if (isTruthy(input$current.cluster) && isTruthy(input$comparison.cluster)) {
    cluster.scatter.plot <- scatter.plot(current.cluster(), comparison.cluster(), cluster.metacells.selected(), cluster.markers(), cluster.markers.selected(), 'cluster')
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.cluster()$exp.label, 
                        current.cluster()$cluster, comparison.cluster()$cluster, cluster.markers()$gene, 
                        cluster.markers.selected()$gene, input$expr.filter.opt, input$pval.thresh,
                        input$max.amt.without, input$min.amt.within))
  } else {
    cluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison cluster in the 'Cluster' panel")    
    key.str <- 'missing_clusters'
  }
  
  height <- img.size.round(session$clientData[[glue("output_gene.expr.scatter.cluster_height")]])
  width <- height
  
  renderCacheImage(cluster.scatter.plot, glue("scatter_cluster_{key.str}"), width, height, progress=progress)
}, deleteFile=FALSE)

output$gene.expr.scatter.cluster.dl <- downloadHandler(filename="scatter.zip", 
                                                       content= function(file) {
                                                         req(isTruthy(input$current.cluster) && isTruthy(input$comparison.cluster))
                                                         cluster.scatter.plot <- scatter.plot(current.cluster(), comparison.cluster(), cluster.metacells.selected(), cluster.markers(), cluster.markers.selected(), 'cluster', return.closure = TRUE)()
                                                         send.zip(cluster.scatter.plot, 'scatter', file, others='scatter-plot.R')
                                                       })

output$gene.expr.scatter.subcluster <- renderImage({
  progress <- shiny.progress('Scatter')
  if (!is.null(progress)) on.exit(progress$close())
  
  if (isTruthy(input$current.subcluster) && isTruthy(input$comparison.subcluster)) {
    subcluster.scatter.plot <- scatter.plot(current.subcluster(), comparison.subcluster(), subcluster.metacells.selected(), subcluster.markers(), subcluster.markers.selected(), 'subcluster')
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.subcluster()$exp.label,
                        current.subcluster()$subcluster, comparison.subcluster()$subcluster, subcluster.markers()$gene,
                        subcluster.markers.selected()$gene, input$expr.filter.opt, input$pval.thresh,
                        input$max.amt.without, input$min.amt.within))
  } else {
    subcluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison subcluster in the 'Cluster' panel")    
    key.str <- 'missing_subclusters'
  }
  
  height <- img.size.round(session$clientData[[glue("output_gene.expr.scatter.subcluster_height")]])
  width <- height
  
  renderCacheImage(subcluster.scatter.plot, glue("scatter_subcluster_{key.str}"), width, height, progress=progress)
}, deleteFile=FALSE)

output$gene.expr.scatter.subcluster.dl <- downloadHandler(filename="scatter.zip", 
                                                          content= function(file) {
                                                            req(isTruthy(input$current.subcluster) && isTruthy(input$comparison.subcluster))
                                                            subcluster.scatter.plot <- scatter.plot(current.subcluster(), comparison.subcluster(), subcluster.metacells.selected(), subcluster.markers(), subcluster.markers.selected(), 'subcluster', return.closure = TRUE)()
                                                            send.zip(subcluster.scatter.plot, 'scatter', file, others='scatter-plot.R')
                                                          })

# Creates a dot plot showing the relative normalized expression for the top.N clusters or subclusters.
# Doesn't work well for multiple genes, because the topN are different and the order also varies.
rank.plot <- function(clusters, kind, genes, return.closure=FALSE) {
  require(tidyr)
  
  Kind <- paste0(toupper(substring(kind,1,1)),substring(kind,2))
  clusters <- gather(clusters, gene.var, value, contains('_target.sum'))
  clusters <- separate(clusters, gene.var, c('gene','var'), sep='_')
  clusters <- spread(clusters, var, value)
  clusters$cx.disp <- (
    if (input$opt.rank.by.region) {
      clusters[[paste0(kind,'.disp')]]
    } else {
      glue("{clusters[[paste0(kind,'.disp')]]} ({clusters$region.abbrev})")
    }
  )
  clusters$gene <- factor(clusters$gene, levels=genes)
  
  clusters <- arrange(filter(clusters, target.sum.per.100k > 0), desc(target.sum.per.100k))
  
  if (nrow(clusters) == 0) { return(plot.text("No data")) }
  
  clusters$cx.disp <- with(clusters, factor(cx.disp, levels=rev(unique(cx.disp))))
  
  clusters.top <- (
    if (input$opt.rank.by.region) {
      group_by(clusters, region.disp,gene) %>% top_n(as.integer(input$top.N), target.sum.per.100k)       
    } else {
      group_by(clusters, gene) %>% top_n(as.integer(input$top.N), target.sum.per.100k)      
    }
  )
  gene.description <- paste(clusters.top$gene,"-",gene.desc(clusters.top$gene))
  clusters.top$gene.description <- factor(gene.description, levels=unique(gene.description[order(clusters.top$gene)]))
  
  rank.facet_grid <- (
    if (input$opt.rank.by.region) {
      facet_grid(region.disp~gene.description, scales="free_y")
    } else {
      facet_grid(~gene.description, scales="free_y")
    }
  )
  
  p.func <- function() {
    source("rank-plot.R", local = TRUE)
    plot.gg   # Not clear why this is required here, but not in other *-plot. 
  }
  
  if (return.closure) {
    p.func
  } else {
    p.func()
  }
}


heatmap.clusters <- reactive({
  select(clusters.selected(), 'region.abbrev','cluster.disp','c.id', contains('_target.sum.per.100k'), -class.disp) %>% unique
})
heatmap.subclusters <- reactive({
  select(subclusters.selected(), 'region.abbrev','subcluster.disp','sc.id', contains('_target.sum.per.100k'), -class.disp) %>% unique
})

observeEvent(input$gene.expr.heatmap.cluster_rows_selected, {
  updateTabsetPanel(session, "controltabs", selected="controltabs-compare")
  updateSelectInput(session, "current.cluster", selected=heatmap.clusters()$c.id[input$gene.expr.heatmap.cluster_rows_selected])  
})

observeEvent(input$gene.expr.heatmap.subcluster_rows_selected, {
  updateTabsetPanel(session, "controltabs", selected="controltabs-compare")
  updateSelectInput(session, "current.subcluster", selected=heatmap.subclusters()$sc.id[input$gene.expr.heatmap.subcluster_rows_selected])  
})

expr.heatmap <- function(cx, kind, genes=user.genes(), opt.heatmap.max=input$opt.heatmap.max.per100k) {
  Kind <- paste0(toupper(substring(kind,1,1)),substring(kind,2))

  gene.cols.idx <- 4:ncol(cx)
  coldefs <- list(list(targets=0:1, 
                       width='40px',
                       render = JS(
                         "function(data, type, row, meta) {",
                         "return type === 'display' && data.length > 40 ?",
                         "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
                         "}")),
                  list(targets=2, visible=FALSE),
                  list(targets=gene.cols.idx-1,
                       class = "dt-center",
                       render = JS("function(data, type, row, meta) {",
                                   "return type==='sort' ? -Math.round(data) :",
                                   "(type==='display' ? '<span title=\"'+Math.round(data)+'\"><img width=\"10px\" height=\"10px\" class=\"img-circle\" style=\"opacity:'+(Math.min(data,",opt.heatmap.max,")/",opt.heatmap.max,")+'\" src=\"data:image/gif;base64,R0lGODlhAQABAIAAAHd3dwAAACH5BAAAAAAALAAAAAABAAEAAAICRAEAOw==\"></span>':data);}"))
                  )
  DT::datatable(cx,
                selection = "single",
                rownames = FALSE,
                colnames = c('Region',Kind, '', genes),
                class = c('cell-border compact'),
                options = list(dom="t", paging=FALSE, autoWidth=TRUE, columnDefs=coldefs))
}

output$gene.expr.heatmap.cluster <- DT::renderDataTable({
  expr.heatmap(heatmap.clusters(),'cluster')
})
output$gene.expr.heatmap.cluster.dl <- downloadHandler(filename="cluster_heatmap.csv",
                                                       content=function(file) { write.csv(heatmap.clusters(), file=file, row.names=FALSE) })

output$gene.expr.heatmap.subcluster <- DT::renderDataTable({
  expr.heatmap(heatmap.subclusters(),'subcluster')
})
output$gene.expr.heatmap.subcluster.dl <- downloadHandler(file="subcluster_heatmap.csv",
                                                          content=function(file) { write.csv(heatmap.subclusters(), file=file, row.names=FALSE) })

opt.plot.height <- reactive({
  if (input$opt.plot.height=='fixed') {
    500
  } else {
    region.mult <- (if (input$opt.rank.by.region) nrow(regions.selected()) else 1)
    gene.mult <- length(user.genes())
    max(500, 16 * region.mult * gene.mult * as.integer(input$top.N))
  }
})

rank.imageOutput <- function(nm) {
  span(class="img-center",imageOutput(nm, height=sprintf("%dpx",opt.plot.height())))
}

output$gene.expr.rank.cluster.output <- renderUI({
  rank.imageOutput("gene.expr.rank.cluster")
})

output$gene.expr.rank.subcluster.output <- renderUI({
  rank.imageOutput("gene.expr.rank.subcluster")
})

output$gene.expr.rank.cluster <- renderPlot({
  if (isTruthy(user.genes())) {
    if (length(user.genes()) > 2) {
      plot.text("Rank display requires 1 or 2 genes.")
    } else {
      clusters <- select(clusters.selected(), -class.disp) %>% unique
      rank.plot(clusters,'cluster', user.genes())
    }
  } else {
    plot.text("Enter a gene symbol in the 'Query' panel\nto display a ranked order of clusters by transcript abundance")    
  }
})

output$gene.expr.rank.subcluster <- renderPlot({
  if (isTruthy(user.genes())) {
    if (length(user.genes()) > 2) {
      plot.text("Rank display requires 1 or 2 genes.")
    } else {
      subclusters <- select(subclusters.selected(), -class.disp) %>% unique
      rank.plot(subclusters, 'subcluster', user.genes())
    }
  } else {
    plot.text("Enter a gene symbol in the 'Query' panel\nto display a ranked order of subclusters by transcript abundance")    
  }
})

output$gene.expr.rank.cluster.dl <- downloadHandler(filename="rank.zip", 
                                                    content= function(file) {
                                                      req(isTruthy(user.genes()))
                                                      
                                                      rank.plot.func <- rank.plot(select(clusters.selected(), -class.disp) %>% unique, 'cluster', user.genes(), return.closure=TRUE)
                                                      send.zip(rank.plot.func, 'rank', file, others='rank-plot.R')
                                                    })

output$gene.expr.rank.subcluster.dl <- downloadHandler(filename="rank.zip", 
                                                       content= function(file) {
                                                         req(isTruthy(user.genes()))
                                                         rank.plot.func <- rank.plot(select(subclusters.selected(), -class.disp) %>% unique, 'subcluster', user.genes(), return.closure=TRUE)
                                                         send.zip(rank.plot.func, 'rank', file, others='rank-plot.R')
                                                       })
