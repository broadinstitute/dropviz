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
# If user.genes exist and option to display cells in gene search, but there is no comparison, yet, then return a stub with user.genes
cluster.markers.selected <- reactive({
  log.reactive("fn: cluster.markers.selected")
  if (isTruthy(input$dt.cluster.markers_rows_selected)) {
    cluster.markers()[input$dt.cluster.markers_rows_selected,]
  } else {
    tibble(gene=character(0))
  }
})

subcluster.markers.selected <- reactive({
  log.reactive("fn: subcluster.markers.selected")
  if (isTruthy(input$dt.subcluster.markers_rows_selected)) {
    subcluster.markers()[input$dt.subcluster.markers_rows_selected,]
  } else {
    tibble(gene=character(0))
  }
})

#####################################################################################################
# OUTPUT - Markers tables

dt.markers <- function(mrkrs,label) {
  if (nrow(mrkrs)==0) {
    rows.selected <- NULL
  } else {
    prev.rows.selected <- isolate(input[[paste0(label,'_rows_selected')]])
    if (!exists(label) || identical(get(label),mrkrs)) {
      rows.selected <- prev.rows.selected
    } else {
      rows.selected <- which(mrkrs$gene %in% get(label)$gene[prev.rows.selected])
    }
  }
  assign(label, mrkrs, 1) # assign global
  DT::datatable(select(mrkrs, gene, description, log.target.u, log.comparison.u, fc.disp, pval, row.highlight),
                rownames = FALSE,
                selection=list(mode="multiple", selected=rows.selected),
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
  dt.markers(mrkrs,'dt.cluster.markers')
})


output$dt.cluster.markers.dl <- downloadHandler(filename="cluster-markers.csv", 
                                                content= function(file) {
                                                  write.csv(cluster.markers(), file=file, row.names=FALSE)
                                                })

output$dt.cluster.markers.heading <- renderUI({
  if (isTruthy(current.cluster.i())) {
    target.names.tbl <- select(inner_join(experiments, current.cluster(), by='exp.label'), exp.abbrev, cluster.disp)
    target.names <- paste(glue("{target.names.tbl$exp.abbrev} {target.names.tbl$cluster.disp}"),collapse='+')
    if (nrow(comparison.cluster())==1 && comparison.cluster()$cluster=='global') {
      comparison.names <- comparison.cluster()$cluster.disp
    } else {
      comparison.names.tbl <- select(inner_join(experiments, comparison.cluster(), by='exp.label'), exp.abbrev, cluster.disp)
      comparison.names <- paste(glue("{comparison.names.tbl$exp.abbrev} {comparison.names.tbl$cluster.disp}"),collapse='+')
    }
    tags$h4(glue("Differentially Over-Expressed: {target.names} vs {comparison.names}"))
  } else {
    tags$p(align="center","Choose a target and comparison cluster in the 'Cluster' panel to find differentially expressed genes")
  }
})

output$dt.subcluster.markers <- DT::renderDataTable( {
  mrkrs <- mutate(subcluster.markers(), row.highlight=ifelse(user.selected,ifelse(expr.pass,1,0),2), description=gene.desc(gene))
  dt.markers(mrkrs,'dt.subcluster.markers')
})


output$dt.subcluster.markers.dl <- downloadHandler(filename="subcluster-markers.csv", 
                                                   content= function(file) {
                                                     write.csv(subcluster.markers(), file=file, row.names=FALSE)
                                                   })

output$dt.subcluster.markers.heading <- renderUI({
  if (isTruthy(current.subcluster.i())) {
    target.names.tbl <- select(inner_join(experiments, current.subcluster(), by='exp.label'), exp.abbrev, subcluster.disp)
    target.names <- paste(glue("{target.names.tbl$exp.abbrev} {target.names.tbl$subcluster.disp}"),collapse='+')
    if (nrow(comparison.subcluster())==1 && comparison.subcluster()$subcluster=='global') {
      comparison.names <- comparison.subcluster()$subcluster.disp
    } else {
      comparison.names.tbl <- select(inner_join(experiments, comparison.subcluster(), by='exp.label'), exp.abbrev, subcluster.disp)
      comparison.names <- paste(glue("{comparison.names.tbl$exp.abbrev} {comparison.names.tbl$subcluster.disp}"),collapse='+')
    }
    
    
    tags$h4(glue("Differentially Over-Expressed: {target.names} vs {comparison.names}"))
  } else {
    tags$p(align="center","Choose a target and comparison subcluster in the 'Cluster' panel to find differentially expressed genes")
  }
})

