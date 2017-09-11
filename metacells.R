source("metacells-shared.R", local=TRUE)

# adds local attributes to results of compute.pair
cx.pairwise <- function(exp.label, cx, cmp.cx, kind) {

  mutate(compute.pair(exp.label, cx, cmp.cx, kind), 
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
  
  cx <- as.character(current.cluster()$cluster)
  cmp.cx <- ifelse(comparison.cluster()$cluster == 'global', paste0('N', current.cluster()$cluster), as.character(comparison.cluster()$cluster))

  cx.pairwise(current.cluster()$exp.label, cx, cmp.cx, 'cluster') %>% 
    mutate(region.disp=current.cluster()$region.disp)  
})

# same as above, but for subclusters
subcluster.metacells.selected <- reactive({
log.reactive("fn: subcluster.metacells.selected")
  req(input$current.subcluster)

  cx <- as.character(current.subcluster()$subcluster)
  cmp.cx <- ifelse(comparison.subcluster()$subcluster == 'global', paste0('N', current.subcluster()$subcluster), as.character(comparison.subcluster()$subcluster))

  cx.pairwise(current.subcluster()$exp.label, cx, cmp.cx, 'subcluster') %>% 
    mutate(region.disp=current.subcluster()$region.disp)
})

# returns a closure of a scatter plot of target vs comparison
scatter.plot <- function(target, comparison, cx.metacells, cx.markers, cx.selected.markers, kind, return.closure=FALSE) {
  function(progress=NULL) {
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
    cluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison cluster in the 'Filter Cells' panel")    
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
    subcluster.scatter.plot <- function(progress) plot.text("Choose a target and comparison subcluster in the 'Filter Cells' panel")    
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
rank.plot <- function(clusters, kind) {
  # TODO: clean up this ugliness
  require(tidyr)
  Kind <- paste0(toupper(substring(kind,1,1)),substring(kind,2))
  clusters <- gather(clusters, gene, amount, ends_with('.log.target.u'))
  clusters$gene <- sub('(.*).log.target.u',"\\1", clusters$gene)
  cx.disp <- paste0(kind,'.disp')
  clusters$cx.disp <- paste(clusters$region.disp,clusters[[cx.disp]])
  clusters <- arrange(clusters, desc(amount))
  clusters$cx.disp <- with(clusters, factor(cx.disp, levels=rev(unique(cx.disp))))
  clusters.top <- group_by(clusters, region.disp,gene) %>% top_n(as.integer(input$top.N), amount) %>%
    mutate(gene.description=paste(gene,"-",gene.desc(gene)))
  ggplot(clusters.top, aes(x=amount, y=cx.disp)) + geom_point(size=3) + 
    ggtitle(glue("Ranked {Kind} by Gene Expression")) + 
    xlab("Normalized log mean") + ylab("") + facet_grid(region.disp~gene.description, scales = "free_y") + coord_cartesian(xlim=c(0,6))
  
}

output$gene.expr.rank.cluster <- renderPlot({
  if (isTruthy(user.genes())) {
    clusters <- select(clusters.selected(), -class.disp) %>% unique
    rank.plot(clusters,'cluster')
  } else {
    plot.text("Enter a gene symbol in the 'Filter Cells' panel\nto display a ranked order of clusters by transcript abundance")    
  }
})

output$gene.expr.rank.subcluster <- renderPlot({
  if (isTruthy(user.genes())) {
    subclusters <- select(subclusters.selected(), -class.disp) %>% unique
    rank.plot(subclusters, 'subcluster')
  } else {
    plot.text("Enter a gene symbol in the 'Filter Cells' panel\nto display a ranked order of subclusters by transcript abundance")    
  }
})
