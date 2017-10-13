# Manages the choice of the current cluster and subcluster for analysis

# There are two parts - 
#  1. server side tracking of the selection ([sub]current.cluster.row)
#  2. output UI elements current.[sub]cluster and comparison.[sub]cluster 


######################################################################
# Current cluster and subcluster (could be factored)

# the index into the selection
current.cluster.i <- reactive({
log.reactive("fn: current.cluster.i")
  as.numeric(input$current.cluster)
})

# the row corresponding to the selection
current.cluster <- reactive({
log.reactive("fn: current.cluster")
  clusters.selected()[current.cluster.i(),]
})

output$current.cluster <- renderUI({
  
  choices <- setNames(1:nrow(clusters.selected()), glue("{clusters.selected()$region.abbrev} {clusters.selected()$cluster.disp}"))
  
  if (length(choices)==1) {
    selected <- choices
  } else {
    selected <- NULL
  }
  
  column(12, selectizeInput("current.cluster", "Target cluster", choices=c("Select target cluster"="",choices), selected=selected, multiple=input$compare.multiple))
})

output$comparison.cluster <- renderUI({
  if (isTruthy(input$current.cluster)) {
    choices <- setdiff(1:nrow(clusters.selected()), current.cluster.i())
    names(choices) <- glue("{clusters.selected()$region.abbrev[choices]} {clusters.selected()$cluster.disp[choices]}")
    
    ## prepend the region comparison
    if (!input$compare.multiple) {
      choices <- c(0,choices) 
      names(choices)[1] <- glue("Rest of {current.cluster()$region.disp}")
    }
    
    column(12, selectizeInput("comparison.cluster", "Comparison", choices=choices, multiple=input$compare.multiple))
  }
})


current.subcluster.i <- reactive({
  log.reactive("fn: current.subcluster.i")
  as.numeric(input$current.subcluster)
})

current.subcluster <- reactive({
  log.reactive("fn: current.subcluster") 
  subc.s <- subclusters.selected_() %>% gene.cols('subcluster') # HACK: Fix Me
  subc.s[current.subcluster.i(),]  
})

output$current.subcluster <- renderUI({

  subc.s <- subclusters.selected_() %>% gene.cols('subcluster') # HACK: Fix Me
  choices <- setNames(1:nrow(subc.s), glue("{subc.s$region.abbrev} {subc.s$subcluster.disp}"))
  
  if (length(choices)==1) {
    selected <- choices
  } else {
    selected <- NULL
  }
  
  column(12, selectizeInput("current.subcluster", "Target Subcluster", choices=c("Select target subcluster"="",choices), selected=selected, multiple=input$compare.multiple))
})

output$comparison.subcluster <- renderUI({
  if (isTruthy(input$current.subcluster)) {
    
    choices <- setdiff(1:nrow(subclusters.selected()), current.subcluster.i())
    names(choices) <- glue("{subclusters.selected()$region.abbrev[choices]} {subclusters.selected()$subcluster.disp[choices]}")
    
    # prepend the region comparison
    if (!input$compare.multiple) {
      choices <- c(0,choices) 
      names(choices)[1] <- glue("Rest of {current.subcluster()$region.disp}")
    }
    
    column(12, selectizeInput("comparison.subcluster", "Comparison", choices=choices, multiple=input$compare.multiple))
  }
})



comparison.cluster <- reactive({
log.reactive("fn: comparison.cluster")
  req(input$comparison.cluster)
  if (length(input$comparison.cluster)==1 && input$comparison.cluster==0) {
    tibble(cluster='global', exp.label=first(current.cluster()$exp.label), cluster.disp=glue("Rest of {current.cluster()$region.disp}"))
  } else {
    clusters.selected()[as.numeric(input$comparison.cluster),]  
  }
})

comparison.subcluster <- reactive({
log.reactive("fn: comparison.subcluster")
  req(input$comparison.subcluster)
  if (length(input$comparison.subcluster)==1 && input$comparison.subcluster==0) {
    tibble(subcluster='global', exp.label=first(current.subcluster()$exp.label), subcluster.disp=glue("Rest of {current.subcluster()$region.disp}"))
  } else {
    subclusters.selected()[as.numeric(input$comparison.subcluster),]  
  }
})

#####################
# Convenience buttons to cycle to next or previous subcluster or cluster
# If the user changes the pull-down, then update the table selection accordingly
observeEvent(input$current.cluster,{
  req(input$current.cluster)
  dt.clusters.proxy %>% DT::selectRows(current.cluster.i())
})

# If the user changes the pull-down, then update the table selection accordingly
observeEvent(input$current.subcluster,{
  req(input$current.subcluster)
  dt.subclusters.proxy %>% DT::selectRows(current.subcluster.i())
})

