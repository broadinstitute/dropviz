# Manages the choice of the current cluster and subcluster for analysis

# There are two parts - 
#  1. server side tracking of the selection ([sub]current.cluster.row)
#  2. output UI elements current.[sub]cluster and comparison.[sub]cluster 


######################################################################

# the index into the selection
current.cluster.i <- reactive({
  if (is.null(input$current.cluster)) NA_integer_ else as.numeric(input$current.cluster)
})

# the row corresponding to the selection
current.cluster <- reactive({
  filter(clusters.selected_(), c.id %in% current.cluster.i())
})

output$current.cluster <- renderUI({
  choices <- setNames(clusters.selected_()$c.id, glue("{clusters.selected_()$region.abbrev} {clusters.selected_()$cluster.disp}"))

  column(12, selectizeInput("current.cluster", "Target cluster",
                            choices=c("Select target cluster"="",choices),
                            selected=if (length(choices)==1) choices else input$current.cluster,
                            multiple=input$compare.multiple))
})

output$comparison.cluster <- renderUI({
  if (nrow(current.cluster())>0) {
    choices <- setdiff(clusters.selected_()$c.id, current.cluster.i())
    choices.nms <- filter(unique(select(clusters.selected_(), c.id, region.abbrev, cluster.disp)), c.id %in% choices)
    names(choices) <- glue("{choices.nms$region.abbrev} {choices.nms$cluster.disp}")
    
    ## prepend the region comparison
    if (!input$compare.multiple) {
      choices <- c(0,choices) 
      names(choices)[1] <- glue("Rest of {current.cluster()$region.disp}")
    }
    
    column(12, selectizeInput("comparison.cluster", "Comparison", 
                              choices=choices,
                              selected=input$comparison.cluster,
                              multiple=input$compare.multiple))
  }
})


current.subcluster.i <- reactive({
  if (is.null(input$current.subcluster)) NA_integer_ else as.numeric(input$current.subcluster)
})

current.subcluster <- reactive({
  subc.s <- subclusters.selected_()
  filter(subc.s, sc.id %in% current.subcluster.i())
})

output$current.subcluster <- renderUI({

  subc.s <- subclusters.selected_()
  choices <- setNames(subc.s$sc.id, glue("{subc.s$region.abbrev} {subc.s$subcluster.disp}"))
  
  column(12, selectizeInput("current.subcluster", "Target Subcluster",
                            choices=c("Select target subcluster"="",choices),
                            selected=if (length(choices)==1) choices else input$current.subcluster,
                            multiple=input$compare.multiple))
})

output$comparison.subcluster <- renderUI({
  if (nrow(current.subcluster())>0) {
    
    choices <- setdiff(subclusters.selected_()$sc.id, current.subcluster.i())
    choices.nms <- filter(unique(select(subclusters.selected_(), sc.id, region.abbrev, subcluster.disp)), sc.id %in% choices)
    names(choices) <- glue("{choices.nms$region.abbrev} {choices.nms$subcluster.disp}")
    
    # prepend the region comparison
    if (!input$compare.multiple) {
      choices <- c(0,choices) 
      names(choices)[1] <- glue("Rest of {current.subcluster()$region.disp}")
    }
    
    column(12, selectizeInput("comparison.subcluster", "Comparison",
                              choices=choices,
                              selected=input$comparison.subcluster,
                              multiple=input$compare.multiple))
  }
})



comparison.cluster <- reactive({
  req(input$comparison.cluster)
  if (length(input$comparison.cluster)==1 && input$comparison.cluster==0) {
    tibble(cluster='global', exp.label=first(current.cluster()$exp.label), cluster.disp=glue("Rest of {current.cluster()$region.disp}"))
  } else {
    filter(clusters.selected_(), c.id %in% as.numeric(input$comparison.cluster))
  }
})

comparison.subcluster <- reactive({
  req(input$comparison.subcluster)
  if (length(input$comparison.subcluster)==1 && input$comparison.subcluster==0) {
    tibble(subcluster='global', exp.label=first(current.subcluster()$exp.label), subcluster.disp=glue("Rest of {current.subcluster()$region.disp}"))
  } else {
    filter(subclusters.selected_(), sc.id %in% as.numeric(input$comparison.subcluster))  
  }
})

#####################
# If the user changes the pull-down, then update the table selection accordingly
observeEvent(input$current.cluster,{
  req(input$current.cluster)
  # FIXME: current.cluster.i is now the ID
#  dt.clusters.proxy %>% DT::selectRows(current.cluster.i())
})

# If the user changes the pull-down, then update the table selection accordingly
observeEvent(input$current.subcluster,{
  req(input$current.subcluster)
  # FIXME: current.cluster.i is now the ID
#  dt.subclusters.proxy %>% DT::selectRows(current.subcluster.i())
})

