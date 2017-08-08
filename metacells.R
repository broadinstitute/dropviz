# return a tibble of cumulative gene expression for the cells in the current.cluster() and comparison.cluster() 
# if the comparison is global, then expression is all cells except the target in the target region.
cluster.metacells.selected <- reactive({
log.reactive("fn: cluster.metacells.selected")
  req(input$current.cluster)
  ddply(current.cluster(), .(exp.label), function(r) {
    fn <- glue("www/metacells/{r$exp.label}.cluster.metacells.RDS")
    readRDS(fn)
  }) %>% as_tibble
})

subcluster.metacells.selected <- reactive({
log.reactive("fn: subcluster.metacells.selected")
  req(input$current.subcluster)
  ddply(current.subcluster(), .(exp.label), function(r) {
    fn <- glue("www/metacells/{r$exp.label}.subcluster.metacells.RDS")
    readRDS(fn)
  }) %>% as_tibble
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
    
    cell.pairs <- cx.metacells[,c('exp.label','gene',as.character(target.cx),as.character(compare.cx))] %>%
      inner_join(region.names(), by='exp.label')
    colnames(cell.pairs) <- c('Region','Gene','Target','Comparison','region.disp')
    cell.pairs <- na.omit(filter(cell.pairs, Target!=0 & Comparison != 0))

    abline <- lm(log10(Comparison) ~ log10(Target), cell.pairs)
    abline.intercept <- coef(abline)[1]
    abline.slope <- coef(abline)[2]
    
    cell.pairs$fc <- (log10(cell.pairs$Comparison) - abline.intercept) / log10(cell.pairs$Target)
    cell.pairs$fc.thresh <- cell.pairs$fc > abline.slope*input$fold.change | cell.pairs$fc < 1/(abline.slope*input$fold.change)
    cell.pairs$filter.thresh <- cell.pairs$Gene %in% cx.markers$GENE
    cell.pairs$user.gene <- cell.pairs$Gene %in% input$user.genes
    cell.pairs$gene.selected <- cell.pairs$Gene %in% cx.selected.markers$GENE
    cell.pairs$color <- factor(with(cell.pairs, ifelse(user.gene, 'User Set', ifelse(gene.selected, 'Selected', ifelse(filter.thresh, 'Pass Criteria', ifelse(fc.thresh, 'Pass FC', 'Other'))))), levels=c('Other','Pass FC','Pass Criteria','Selected','User Set'))
    cell.pairs$size <- factor(ifelse(cell.pairs$color %in% c('Selected','User Set','Pass Criteria'),'Filtered','Other'), levels=c('Other','Filtered'))
    cell.pairs$gene.label <- ifelse(cell.pairs$gene.selected | cell.pairs$user.gene, cell.pairs$Gene, NA_character_)
    
    opt.scatter.gene.labels <- input$opt.scatter.gene.labels
    opt.expr.size <- input$opt.expr.size
    p.func <- function() {
      require(ggplot2); require(dplyr)
      if (opt.scatter.gene.labels && any(!is.na(cell.pairs$gene.label))) {
        gene.labels <- geom_label(aes(x=Target*1.5, label=gene.label), size=3, show.legend = FALSE)
      } else {
        gene.labels <- geom_blank()
      }
      
      plot.gg <- ggplot(dplyr::filter(cell.pairs, Target>0 & Comparison>0), aes(x=Target, y=Comparison, color=color, size=size)) + geom_point() + 
        geom_abline(intercept=abline.intercept, slope=abline.slope, color='pink') + gene.labels +
        scale_x_log10() + scale_y_log10() + scale_color_manual(guide="none",values=c('Other'='grey','Pass FC'='lightgreen','Pass Criteria'='blue','Selected'='purple','User Set'='black')) +
        scale_size_manual(guide="none",values=c('Other'=1,'Filtered'=opt.expr.size)) + facet_wrap(~region.disp) + xlab(target.label) + ylab(compare.label)

      if (interactive()) {
        print(plot.gg)
      } 
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
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.cluster()$exp.label, current.cluster()$cluster, comparison.cluster()$cluster, cluster.markers()$GENE, cluster.markers.selected()$GENE))
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
    key.str <- digest(c(input$fold.change, input$opt.scatter.gene.labels,current.subcluster()$exp.label, current.subcluster()$subcluster, comparison.subcluster()$subcluster, subcluster.markers()$GENE, subcluster.markers.selected()$GENE))
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
