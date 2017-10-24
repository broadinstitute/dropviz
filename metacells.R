source("metacells-shared.R", local=TRUE)

# in prep-metacells, all cx vs Ncx are pre-computed. Copy those locally. This does not work on Windows
try(system2("rsync",c("--size-only",glue("{prep.dir}/pairs/*"),glue("{cache.dir}/metacells/"))))

# adds local attributes to results of compute.pair
cx.pairwise <- function(exp.label, cx, cmp.exp.label, cmp.cx, kind) {

  mutate(compute.pair(exp.label, cx, cmp.exp.label, cmp.cx, kind), 
         fc.thresh=( if (input$expr.filter.opt %in% 'fc') { (abs(fc) > log(input$fold.change)) } else { fc > log(input$fold.change) } ),
         pval.thresh=pval < 10^input$pval.thresh,
         amt.thresh=( (log.target.u > input$min.amt.within) & (log.comparison.u < input$max.amt.without) ),
         user.selected=gene %in% user.genes(),
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
  req(current.cluster()$cluster)
  
  cx <- as.character(current.cluster()$cluster)
  cmp.cx <- ifelse(comparison.cluster()$cluster == 'global', paste0('N', current.cluster()$cluster), as.character(comparison.cluster()$cluster))

  cx.pairwise(current.cluster()$exp.label, cx, comparison.cluster()$exp.label, cmp.cx, 'cluster')
})

# same as above, but for subclusters
subcluster.metacells.selected <- reactive({
log.reactive("fn: subcluster.metacells.selected")
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
                         target.region=paste(target$region.disp, collapse=' + '),
                         compare.region=paste(comparison$region.disp, collapse=' + ')
    )


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
          fc.line1 <- geom_abline(intercept=log(input$fold.change), slope=1, alpha=0.5, color='grey') 
        } else {
          fc.line1 <- geom_blank()
        }
        fc.line2 <- geom_abline(intercept=-log(input$fold.change), slope=1, alpha=0.5, color='grey')
      } else {
        fc.line1 <- fc.line2 <- geom_blank()
      }
      
      if (input$expr.filter.opt %in% c('both','amt','either')) {
        amt.line1 <- geom_hline(yintercept = input$max.amt.without, color='grey') 
        amt.line2 <- geom_vline(xintercept = input$min.amt.within, color='grey')
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
        facet_grid(compare.region~target.region) + xlab(paste(target.label,"\n(log normal mean)")) + ylab(paste(compare.label,"\n(log normal mean)"))

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
    cluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison cluster in the 'Compare' panel")    
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
                                                         send.zip(cluster.scatter.plot, 'scatter', file)
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
    subcluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison subcluster in the 'Compare' panel")    
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
                                                            send.zip(subcluster.scatter.plot, 'scatter', file)
                                                          })

# Creates a dot plot showing the relative normalized expression for the top.N clusters or subclusters.
# Doesn't work well for multiple genes, because the topN are different and the order also varies.
rank.plot <- function(clusters, kind, genes) {
  require(tidyr)
  
  Kind <- paste0(toupper(substring(kind,1,1)),substring(kind,2))
  clusters <- gather(clusters, gene.var, value, contains('-target.sum'))
  clusters <- separate(clusters, gene.var, c('gene','var'), sep='-')
  clusters <- spread(clusters, var, value)
  clusters$cx.disp <- paste(clusters$region.disp,clusters[[paste0(kind,'.disp')]])
  clusters$gene <- factor(clusters$gene, levels=genes)
  
  clusters <- arrange(clusters, desc(target.sum.per.100k))
  clusters$cx.disp <- with(clusters, factor(cx.disp, levels=rev(unique(cx.disp))))
  
  clusters.top <- group_by(clusters, region.disp,gene) %>% top_n(as.integer(input$top.N), target.sum.per.100k)
  gene.description <- paste(clusters.top$gene,"-",gene.desc(clusters.top$gene))
  clusters.top$gene.description <- factor(gene.description, levels=unique(gene.description[order(clusters.top$gene)]))

  rank.facet_grid <- (
    if (input$opt.rank.by.region) {
      facet_grid(region.disp~gene.description, scales="free_y")
    } else {
      facet_grid(~gene.description, scales="free_y")
    }
  )

  ggplot(clusters.top, aes(x=target.sum.per.100k, xmin=target.sum.L.per.100k, xmax=target.sum.R.per.100k, y=cx.disp, yend=cx.disp)) + geom_point(size=3) + geom_segment(aes(x=target.sum.L.per.100k,xend=target.sum.R.per.100k)) +
    ggtitle(glue("Ranked {Kind}s by Gene Expression")) + 
    xlab(glue("Transcripts Per 100,000 in {Kind}")) + ylab("") + rank.facet_grid 
  
}


heatmap.clusters <- reactive({
  select(clusters.selected(), 'region.abbrev','cluster.disp','c.id', contains('-target.sum.per.100k'), -class.disp) %>% unique
})
heatmap.subclusters <- reactive({
  select(subclusters.selected(), 'region.abbrev','subcluster.disp','sc.id', contains('-target.sum.per.100k'), -class.disp) %>% unique
})

observeEvent(input$gene.expr.heatmap.cluster_rows_selected, {
  updateTabsetPanel(session, "controltabs", selected="controltabs-compare")
  updateSelectInput(session, "current.cluster", selected=heatmap.clusters()$c.id[input$gene.expr.heatmap.cluster_rows_selected])  
})

observeEvent(input$gene.expr.heatmap.subcluster_rows_selected, {
  updateTabsetPanel(session, "controltabs", selected="controltabs-compare")
  updateSelectInput(session, "current.subcluster", selected=heatmap.subclusters()$sc.id[input$gene.expr.heatmap.subcluster_rows_selected])  
})

expr.heatmap <- function(cx, kind, genes=input$user.genes, opt.heatmap.max=input$opt.heatmap.max.per100k) {
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
output$gene.expr.heatmap.subcluster <- DT::renderDataTable({
  expr.heatmap(heatmap.subclusters(),'subcluster')
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
    plot.text("Enter a gene symbol in the 'Compare' panel\nto display a ranked order of clusters by transcript abundance")    
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
    plot.text("Enter a gene symbol in the 'Compare' panel\nto display a ranked order of subclusters by transcript abundance")    
  }
})
