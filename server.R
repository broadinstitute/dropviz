library(digest)
library(ggplot2)
library(shiny)


shiny.progress <- function() {
  tryCatch({
    shiny::Progress$new();
  }, error=function(err) NULL)
}

# uncomment to debug within RStudio
# reactive <- function(expr, env=parent.frame(), ...) exprToFunction(expr, env=env)

shinyServer(function(input, output, session) {
  
  # general plot utility functions
  source("plot.R", local=TRUE)
  
  # display labels for regions, clusters, subclusters
  source("display_labels.R", local=TRUE)
  
  # filtering on [sub]clusters
  source("cell_types.R", local=TRUE)

  # choosing a current [sub]cluster selection among the filtered options
  # for marker filtering and plotting
  source("user_cluster_selection.R", local=TRUE)  
  
  # filtering on differentially expressed genes for the current [sub]cluster selection
  source("markers.R", local=TRUE)    

  # plots for both cluster and subcluster with overlays for labels, gene expression and components
  source("tSNE.R", local=TRUE)  
  
  # differential expression
  source("diffex.R", local=TRUE)
  
  # metacells - scatter plots, heatmaps
  source("metacells.R", local=TRUE)
  
  # independent components
  source("components.R", local=TRUE)

  # proxy objects for shiny
  source("proxy.R", local=TRUE)

  #####################################################################################################
  # save the latest input for interactive debug
  # load with input <- readRDS("dump.RDS")
  observeEvent(input$dump, {
    write.log("Dump"); saveRDS(reactiveValuesToList(input), file="dump.RDS")
  })
  
  observeEvent(input$clear.cache, {
    unlink("www/cache/*.png")
  })

  # dropped this idea of expanding and collapsing the sidebar, but this might be an approach
  # observe({
  #   if (input$wideside1) {
  #     jqui_switch_class("#sidebar", 'col-sm-4', 'col-sm-6',1000)
  #     jqui_switch_class("#mainpanel", 'col-sm-8', 'col-sm-6',1000)
  #   } else {
  #     jqui_switch_class("#sidebar", 'col-sm-6', 'col-sm-4',1000)
  #     jqui_switch_class("#mainpanel", 'col-sm-6', 'col-sm-8',1000)
  #   }
  # })
  
})
