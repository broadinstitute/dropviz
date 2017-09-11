# routines related to independent components

#####################################################################################################
# OUTPUT - Components tables

# returns the major cluster associated with the target subcluster
component.cluster <- reactive({
  log.reactive("fn: component.cluster")
  
  if (isTruthy(input$current.subcluster)) {
    na.omit(current.subcluster())
  } else if (length(unique(subclusters.selected_()$cluster))==1) {
    subclusters.selected_()[1,]
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

# returns a matrix the weights ("rotations") of clusters (columns) on cells (rows) for the clusters.selected.components
component.cell.weights <- reactive({
  log.reactive("fn: component.cell.weights")
  if (nrow(component.cluster())==0)
    return(matrix())
  
  exp <- filter(experiments, exp.label==component.cluster()$exp.label)
  fn <- list.files(glue("{exp$exp.dir}/components"), glue("{exp$base}.cluster{component.cluster()$cluster}\\..*.ICA.RDS"))
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


output$ic.grid <- renderImage({
  fn <- character(0)
  if (nrow(component.cluster())>0) {
    fn <- list.files(glue("{cache.dir}/ic"), glue("{component.cluster()$exp.label}_{component.cluster()$cluster}_{input$opt.components}_[0-9]+_[0-9]+.png"))
  }
  
  if (length(fn)>0) {
    list(src=glue("{cache.dir}/ic/{fn[1]}"))
  } else {
    renderCacheImage(function(progress) plot.text("To display ICs for a specific cluster, choose a target cluster in the 'Filter Cells' panel"),
                     "no-ICs", width = 500, height=500)
  }
}, deleteFile = FALSE)

output$dt.components <- DT::renderDataTable( {
  req(nrow(clusters.selected.components())>0)
  components.tbl <- clusters.selected.components() %>% 
    mutate(Loadings=glue("<img height='75' width='250' src='cache/ic/ic_{exp.label}_{cluster}_IC{ic.number}_250_75.png'/>")) %>%
    dplyr::select(IC=ic.number, Class=cell_class, Status=status, "Annotation Notes"=hypothesized_common_name, Loadings, Region=anatomical_region, Clustering=use_for_clustering)
  
  DT::datatable(components.tbl,
                rownames=FALSE, escape=FALSE,
                selection="multiple",
                options=list(dom="tp", pageLength=50))
})

output$dt.components.heading <- renderUI({
  if (nrow(component.cluster())==1) {
    tags$h4(glue("ICs for {component.cluster()$cluster.disp}"))
  } else {
    tags$p(align="center","To display ICs for a specific cluster, choose a target cluster in the 'Filter Cells' panel")
  }
})

