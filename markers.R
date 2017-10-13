# Differentially expressed genes (i.e. "markers")
#
# Composed of three parts:
# 1. UI elements for specifying expression criteria
# 2. filtering code to retrieve only matching genes
# 3. UI tables for showing gene results

#####################################################################################################
# Most of the UI inputs are static in ui.R. But depending on the expr.filter.opt, the widgets are
# enabled/disabled here.
# 
# The user can choose Fold Change, Amount or BOTH or EITHER
# This disables/enables fold change and amount depending on the user option
observeEvent(input$expr.filter.opt, {
  fc.widgets <- c('fold.change','pval.thresh')
  if (input$expr.filter.opt=='amt') {
    lapply(fc.widgets, function(w) shinyjs::disable(w) )
  } else {
    lapply(fc.widgets, function(w) shinyjs::enable(w) )
  }
  
  amt.present.widgets <- c('min.amt.within','max.amt.without')
  if (input$expr.filter.opt=='fc') {
    lapply(amt.present.widgets, function(w) shinyjs::disable(w) )
  } else {
    lapply(amt.present.widgets, function(w) shinyjs::enable(w) )
  }      
})

#####################################################################################################
# Gene Expression Filter


pairwise.markers <- function(kind) {
  (if (kind=='cluster') {
    cluster.metacells.selected()
  } else {
    subcluster.metacells.selected()
  }) %>% filter(user.selected | expr.pass) %>% arrange(-user.selected, -expr.pass, -fc)
}

# Returns a table of the markers associated with the filtered clusters, 
# filtered according to expression criteria and user genes
cluster.markers <- reactive({
log.reactive("fn: cluster.markers")
  req(current.cluster.i())
  pairwise.markers('cluster')
})


subcluster.markers <- reactive({
log.reactive("fn: subcluster.markers")
  req(current.subcluster.i())
  pairwise.markers('subcluster')
})

# Returns a subset of the cluster.markers for rows that user clicked in table
cluster.markers.selected <- reactive({
log.reactive("fn: cluster.markers.selected")
  if (is.null(input$dt.cluster.markers_rows_selected) || !isTruthy(current.cluster.i())) {
    tibble(gene=character(0))
  } else {
    cluster.markers()[input$dt.cluster.markers_rows_selected,]
  }
})

subcluster.markers.selected <- reactive({
log.reactive("fn: subcluster.markers.selected")
  if (is.null(input$dt.subcluster.markers_rows_selected) || !isTruthy(current.subcluster.i())) {
    tibble(gene=character(0))
  } else {
    subcluster.markers()[input$dt.subcluster.markers_rows_selected,]
  }
})

#####################################################################################################
# OUTPUT - Markers tables

dt.markers <- function(mrkrs) {
  DT::datatable(select(mrkrs, gene, description, log.target.u, log.comparison.u, fc.disp, pval, row.highlight),
                rownames = FALSE,
                selection="multiple",
                colnames = c('Gene','Description', 'Target\n(normalized mean log)', 'Comparison\n(normalized mean log)', 'Fold Ratio', 'P-Value', 'row.highlight'),
                options=list(dom="tp", pageLength=50,
                             language=list(zeroRecords = "No results - adjust Diff Expr criteria using Scatter Plot as guide"),
                             columnDefs = list(list(visible=FALSE, targets=6)))) %>% 
    DT::formatStyle('row.highlight', target='row', 
                    backgroundColor = DT::styleEqual(c(0,1,2), c('pink','lightgreen','white'))) %>%
    DT::formatSignif(c('log.target.u','log.comparison.u','fc.disp','pval'), 3) 
}

output$dt.cluster.markers <- DT::renderDataTable( {
  mrkrs <- mutate(cluster.markers(), row.highlight=ifelse(user.selected,ifelse(expr.pass,1,0),2), description=gene.desc(gene)) 
  dt.markers(mrkrs)
})

output$dt.cluster.markers.dl <- downloadHandler(filename="cluster-markers.csv", 
                                                content= function(file) {
                                                  write.csv(cluster.markers(), file=file)
                                                })

output$dt.cluster.markers.heading <- renderUI({
  if (isTruthy(current.cluster.i())) {
    target.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label%in%current.cluster()$exp.label]} {current.cluster()$cluster.disp}"),collapse='+')
    if (comparison.cluster()$cluster=='global') {
      comparison.names <- comparison.cluster()$cluster.disp
    } else {
      comparison.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label%in%comparison.cluster()$exp.label]} {comparison.cluster()$cluster.disp}"),collapse='+')
    }
    
    
    tags$h4(glue("Differentially Over-Expressed: {target.names} vs {comparison.names}"))
  } else {
    tags$p(align="center","Choose a target and comparison cluster in the 'Compare' panel to find differentially expressed genes")
  }
})

output$dt.subcluster.markers <- DT::renderDataTable( {
  mrkrs <- mutate(subcluster.markers(), row.highlight=ifelse(user.selected,ifelse(expr.pass,1,0),2), description=gene.desc(gene))
  dt.markers(mrkrs)
})

output$dt.subcluster.markers.dl <- downloadHandler(filename="subcluster-markers.csv", 
                                                   content= function(file) {
                                                     write.csv(subcluster.markers(), file=file)
                                                   })

output$dt.subcluster.markers.heading <- renderUI({
  if (isTruthy(current.subcluster.i())) {
    target.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label %in% current.subcluster()$exp.label]} {current.subcluster()$subcluster.disp}"),collapse=' + ')
    if (comparison.subcluster()$subcluster=='global') {
      comparison.names <- comparison.subcluster()$subcluster.disp
    } else {
      comparison.names <- paste(glue("{experiments$exp.abbrev[experiments$exp.label %in% comparison.subcluster()$exp.label]} {comparison.subcluster()$subcluster.disp}"),collapse=' + ')
    }
    
    
    tags$h4(glue("Differentially Over-Expressed: {target.names} vs {comparison.names}"))
  } else {
    tags$p(align="center","Choose a target and comparison subcluster in the 'Compare' panel to find differentially expressed genes")
  }
})

