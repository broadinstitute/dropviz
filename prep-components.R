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
log("Sources loaded")

dir.create(glue("{cache.dir}/ic"), recursive = TRUE, showWarnings=FALSE)


# populate input with the minimum values to generate a faceted tSNE of components

input <- list(opt.components='real', opt.cluster.disp='all', opt.region.disp='region', use.common.name=TRUE)
ic.kinds <- c('real','clustering','all')

no.xy <- theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.title.y=element_blank(),
               axis.text.y=element_blank())

# currently the ICs are generated from the current.subcluster. In the future, it may be parameterized differently.
dlply(experiments, .(exp.label), function(exp) {
  input$tissue <<- exp$exp.title
  ldply(ic.kinds, function(ic.kind) {
    input$opt.components <<- ic.kind
    cluster.idx <- c(1, 1+which(diff(as.numeric(subclusters.selected_()$cluster))!=0))
    cluster.idx <- setNames(cluster.idx, subclusters.selected_()$cluster[cluster.idx])
    ldply(cluster.idx, function(i) {
      input$current.subcluster <<- i
      nr <- nrow(clusters.selected.components())
      if (nr > 0) {
        input$dt.components_rows_selected <<- 1:nr
        sc <- selected.component.cell.weights.xy()
      } else {
        sc <- tibble()
      }
        
      if (nrow(sc)>0) {
        ic.plot <- function() {
          ggplot() + geom_point(data=sc, aes(x=V1, y=V2, color=weight), size=2, alpha=1) + 
            facet_wrap(~region.disp+facet.gg+IC, ncol=4) +
            scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(sc$weight),max(sc$weight)))
        }
        width <- 1000
#        height <- width/2 * ((nr-1)%/%4+1)
        height <- width/2 * (nr %/% 12 + 1)
      } else {
        ic.plot <- function() plot.text("No Data")
        width <- height <- 500
        log(glue("No data for {component.cluster()$cluster}"))
      }
      
      renderCacheImage(ic.plot, glue("ic/ic_{exp$exp.label}_{component.cluster()$cluster}_{ic.kind}"), width, height, opt.use.cache = FALSE)
      
      # create little plots to embed in IC table
      if (ic.kind=='all') {
        ddply(sc, .(IC), function(df) {
          df <- arrange(df, cell)

          red.blue.gradient <-  scale_color_gradient2(guide="none", low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(df$weight),max(df$weight)))
          
          renderCacheImage(function() { ggplot(df, aes(x=cell, y=weight, color=weight)) + geom_point(size=1, alpha=0.7) + no.xy + red.blue.gradient},
                           glue("ic/ic_{exp$exp.label}_{component.cluster()$cluster}_{first(df$IC)}"),
                           width=250, height=75, opt.use.cache = FALSE)
          NULL
        })
      }
      data.frame()
    })
  })
})
