# routines related to independent components

#####################################################################################################
# OUTPUT - Components tables

# returns the major cluster associated with either the target subcluster, or if there is only one
# cluster selected, then use that
component.cluster <- reactive({
  log.reactive("fn: component.cluster")
  
  if (isTruthy(input$current.subcluster)) {
    na.omit(current.subcluster())
  } else {
    tibble()
  }
})

# returns only the components associated with the cluster of the currently selected subcluster
clusters.selected.components <- reactive({
  log.reactive("fn: clusters.selected.components")
  if (nrow(component.cluster())==0)
    return(tibble())
  
  cmp <- filter(components, cluster==component.cluster()$cluster & exp.label==component.cluster()$exp.label)
  
  if (input$opt.components=='real') {
    filter(cmp, status=='Real')
  } else if (input$opt.components=='clustering') {
    filter(cmp, use_for_clustering)
  } else {
    cmp
  }
})

# returns the ICs that the user selected in the table
selected.components <- reactive({
  log.reactive("fn: selected.components")
  if (isTruthy(input$dt.components_rows_selected)) {
    clusters.selected.components()[input$dt.components_rows_selected,]
  } else {
    tibble(ic.number=integer())
  }
})

# returns a matrix the weights ("rotations") on cells for the clusters.selected.components
component.cell.weights <- reactive({
  log.reactive("fn: component.cell.weights")
  if (is.null(component.cluster()))
    return(matrix())
  
  exp <- filter(experiments, exp.label==component.cluster()$exp.label)
  fn <- list.files(glue("{exp$exp.dir}/components"), glue("{exp$base}.cluster{component.cluster()$cluster}.*.ICA.RDS"))
  stopifnot(length(fn)==1)
  readRDS(glue("{exp$exp.dir}/components/{fn}"))$cell_rotations
})

# returns as a tibble the slice of the component.cell.weights matrix for only the components selected in the table
selected.component.cell.weights <- reactive({
  require(tidyr)
  log.reactive("fn: component.cell.weights.selected")
  if (nrow(selected.components())==0)
    return(tibble())
  
  cell.weights <- component.cell.weights()[,selected.components()$ic.number, drop=FALSE]
  colnames(cell.weights) <- paste0("IC", selected.components()$ic.number)  
  bind_cols(tibble(cell=rownames(cell.weights)), as_tibble(cell.weights)) %>% gather(IC, weight, -cell)
})

# adds XY coordinates in local tSNE space to selected.component.cell.weights
selected.component.cell.weights.xy <- reactive({
  if (nrow(selected.components())==0)
    return(tibble())
  
  inner_join(selected.component.cell.weights(), local.xy.selected(), by='cell')
})


#   
#   DT::datatable(select(mrkrs, GENE, FOLDch, pct.1, pct.2, row.highlight),
#                 rownames = FALSE,
#                 selection="multiple",
#                 colnames = c('Gene','Fold Change','Target Present','Other Present','row.highlight'),
#                 options=list(dom="tp", pageLength=50, columnDefs = list(list(visible=FALSE, targets=4)))) %>% 
#     DT::formatStyle('row.highlight', target='row', 
#                     backgroundColor = DT::styleEqual(c(0,1,2), c('pink','lightgreen','white'))) %>%
#     DT::formatSignif('FOLDch', 3) %>% DT::formatPercentage('pct.1', 1) %>% DT::formatPercentage('pct.2', 1)
# }

output$dt.components <- DT::renderDataTable( {
  req(nrow(clusters.selected.components())>0)
  DT::datatable(clusters.selected.components() %>% dplyr::select(IC=ic.number, Class=cell_class, Status=status, Name=hypothesized_common_name, Region=anatomical_region, Clustering=use_for_clustering),
                rownames=FALSE,
                selection="multiple",
                options=list(dom="tp", pageLength=50))
})

