# pre-generate PNGs of the components in tSNE space for each subcluster
# HACK: pretend that we're running shiny and choose inputs so that we cycle through each
# cluster's independent components. This means choosing a subcluster from each cluster and
# then selecting all the IC rows.

library(ggplot2)

# create an environment similar to running a shiny app
output <- list()
source("global.R", local=TRUE)
reactive <- function(expr, env=parent.frame()) exprToFunction(expr, env=env)
source("plot.R", local=TRUE)
source("display_labels.R", local=TRUE)
source("cell_types.R", local=TRUE)
source("user_cluster_selection.R", local=TRUE)  
source("tSNE.R", local=TRUE)  
source("components.R", local=TRUE)
source("markers.R", local=TRUE)
write.log("Sources loaded")

dir.create(glue("{cache.dir}/ic"), recursive = TRUE, showWarnings=FALSE)


# populate input with the minimum values to generate a faceted tSNE of components

input <- list(opt.components='real', opt.cluster.disp='all', opt.region.disp='region', use.common.name=TRUE, opt.downsampling.method="none")
ic.kinds <- c('real','clustering','all')

no.xy <- theme_few() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        strip.text.x=element_text(size=20),
        strip.text.y=element_text(size=14))

# currently the ICs are generated from the current.subcluster. In the future, it may be parameterized differently.
dlply(experiments, .(exp.label), function(exp) {
  input$tissue <<- exp$exp.title
  ldply(ic.kinds, function(ic.kind) {
    input$cell.cluster <<- NULL
    input$opt.components <<- ic.kind
    clusters <- clusters.selected_()
    write.log(glue("{input$tissue} {input$opt.components} {nrow(clusters)} clusters"))
    ddply(clusters, .(c.id), function(cl) {
      input$cell.cluster <<- cl$cluster.disp
      nr <- nrow(clusters.selected.components())
      write.log(nr," components for ", cl$cluster.disp)
      if (nr > 0) {
        input$dt.components_rows_selected <<- 1:nr
        sc.xy <- selected.component.cell.weights.xy()
      } else {
        sc.xy <- tibble()
      }
        
      if (nrow(sc.xy)>0) {
        ic.plot <- function(progress) {
          if (ic.kind == 'all') { kind.text <- 'All' } else if (ic.kind=='real') { kind.text <- '"Real"' } else { kind.text <- '' }
          subcluster.text <- ifelse(ic.kind=='clustering', 'used for subclustering', '')
          title <- glue("{kind.text} ICs for {cl$region.disp} cluster {cl$cluster.disp} {subcluster.text}")
          ggplot() + geom_point(data=sc.xy, aes(x=V1, y=V2, color=weight), size=2, alpha=1) + 
            facet_wrap(~IC, ncol=4) + ggtitle(title) +
            scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(sc.xy$weight),max(sc.xy$weight)))
        }
        width <- 1000
#        height <- width/2 * ((nr-1)%/%4+1)
        height <- width/2 * (nr %/% 12 + 1)
      } else {
        ic.plot <- function(progress) plot.text(glue("No ICs for {cl$region.disp} cluster {cl$cluster.disp}"))
        width <- height <- 500
        write.log(glue("No data for {cl$cluster}"))
      }
      
      renderCacheImage(ic.plot, glue("ic/ic_{exp$exp.label}_{cl$cluster}_{ic.kind}"), width, height, opt.use.cache = FALSE)
      
      # create little plots to embed in IC table
      sc <- selected.component.cell.weights()
      if (ic.kind=='all') {
        ddply(sc, .(IC), function(df) {
          df <- arrange(df, cell)

          red.blue.gradient <-  scale_color_gradient2(guide="none", low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(df$weight),max(df$weight)))
          
          renderCacheImage(function(progress) { ggplot(df, aes(x=cell, y=weight, color=weight)) + geom_point(size=1, alpha=0.7) + no.xy + red.blue.gradient},
                           glue("ic/ic_{exp$exp.label}_{cl$cluster}_{first(df$IC)}"),
                           width=250, height=75, opt.use.cache = FALSE)
          NULL
        })
      }
      data.frame()
    })
  })
})
