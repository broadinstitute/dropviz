library(digest)
library(ggplot2)
library(shiny)

# uncomment to debug within RStudio
#reactive <- function(expr, env=parent.frame()) exprToFunction(expr, env=env)

# On error, return a plot with txt in the center
plot.text <- function(txt) {
  ggplot(data.frame(x=0,y=0,label=txt)) + geom_text(aes(x=x,y=y,label=label)) + theme_void()
}

# adjust to nearest 100. This is to avoid flickering when size of region changes by small amounts in browser
img.size.round <- function(x) {
  (x %/% 100) * 100
}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  # create a plot from plot.func, save it as a PNG using the size of the output region, cache it using key,
  # and return image info, setting the class to the output.id and its id to key.
  renderCacheImage <- function(plot.func, key, width, height=width, opt.use.cache=input$opt.use.cache) {
    # width <- clientData[[glue("output_{output.id}_width")]]
    # height <- clientData[[glue("output_{output.id}_height")]]
    # log(glue("clientData: WxH = {width}x{height}"))
    # width <- 1000
    # height <- 500 * height.hint
    log(glue("WxH = {width}x{height}"))
    
    if (is.null(width) | is.null(height)) {
      stop("Missing width or height}")
    }
    fn <- glue("{cache.dir}/{key}_{width}_{height}.png")
    
    if (!file.exists(fn) || !opt.use.cache) {
      log(glue("Generating plot {fn}"))
      a.plot <- plot.func()
      # log(glue("Printing plot to {fn}"))
      png(fn, width=width, height=height)
      print(a.plot)
      dev.off()
    } else {
      log(glue("Retrieving cached {fn}"))
    }
    
    list(src=fn, width=width, height=height, id=key)
  }
  
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

  #####################################################################################################
  # save the latest input for interactive debug
  # load with input <- readRDS("dump.RDS")
  observeEvent(input$dump, {
    log("Dump"); saveRDS(reactiveValuesToList(input), file="dump.RDS")
  })
  
  observeEvent(input$clear.cache, {
    unlink("www/cache/*.png")
  })
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
