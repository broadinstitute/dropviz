# routines related to independent components

#####################################################################################################
# Components table and tSNE

# essentially just an alias for clusters.selected()
component.cluster <- reactive({

  if (nrow(clusters.selected_()) > 0) {
    clusters.selected_()
  } else {
    tibble()
  }
})

# returns only the components associated with the currently selected cluster(s)
clusters.selected.components <- reactive({
  req(input$opt.components)
  
  if (nrow(component.cluster())==0)
    return(tibble())
  
  cmp <- inner_join(components, component.cluster(), by=c('cluster', 'exp.label'))
  
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
  if (isTruthy(input$dt.components_rows_selected) && nrow(clusters.selected.components())>0) {
   clusters.selected.components()[input$dt.components_rows_selected,]
  } else {
    tibble(ic.number=integer())
  }
})

component.ICA <- reactive({
  if (nrow(component.cluster())==0)
    return(list(cell_rotations=matrix(), genes=matrix()))
  
  validate(need(length(unique(selected.components()$cluster))==1, "Select ICs from only one cluster"))
  
  exp <- filter(experiments, exp.label==first(selected.components()$exp.label))
  fn <- list.files(glue("{exp$exp.dir}/components"), glue("{exp$base}.cluster{selected.components()$cluster}\\..*.ICA.RDS"))
  #  stopifnot(length(fn)==1)
  fn <- fn[1]  # FIXME: multiple matches
  readRDS(glue("{exp$exp.dir}/components/{fn}"))
})

# returns a matrix the weights ("rotations") of clusters (columns) on cells (rows) for the clusters.selected.components
component.cell.weights <- reactive({
  component.ICA()$cell_rotations
})

# returns as a tibble the slice of the component.cell.weights matrix for only the components selected in the table
selected.component.cell.weights <- reactive({
  require(tidyr)
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

# weighs gene transcripts by the cell weight
selected.component.cell.weights.genes <- reactive({
  req(nrow(selected.components()))
  req(!is.null(selected.components()$positive_loading) || !is.null(selected.components()$negative_loading))
  top.genes <- unlist(strsplit(c(selected.components()$positive_loading,selected.components()$negative_loading), '\\s+', perl=TRUE))
  m <- component.ICA()$gene[top.genes, selected.components()$ic.number]
  tibble(gene=names(m), description=gene.desc(names(m)), loading=m) %>% arrange(-loading)
})

output$dt.IC.genes <- DT::renderDataTable({
  DT::datatable(selected.component.cell.weights.genes(),
                rownames=FALSE,
                selection="none",
                options=list(dom="tp")) %>%
    DT::formatSignif('loading')
})

output$dt.IC.subclusters <- DT::renderDataTable({
  req(nrow(selected.component.cell.weights.xy()))
  IC.subc.tbl <- selected.component.cell.weights.xy() %>% group_by(exp.label, cx) %>% 
    summarize(weight=mean(weight)) %>% 
    arrange(-weight) %>% 
    rename("subcluster"="cx", "weight"="weight") %>%
    inner_join(subcluster.names(), by=c('exp.label','subcluster'))
  DT::datatable(IC.subc.tbl[,c("subcluster.disp","weight")], colnames = c("Subcluster","Weight"), rownames=FALSE, selection="none", options=list(dom="tp")) %>%
    DT::formatSignif('weight')
})

output$tsne.selected.component <- renderImage({
  progress <- shiny.progress('IC t-SNE')
  if (!is.null(progress)) on.exit(progress$close())

  sc <- selected.components()
  req(nrow(sc)==1)
  sc.weights <- selected.component.cell.weights.xy()

  tsne.plot <- function(progress) {
    ggplot(sc.weights, aes(x=V1, y=V2, color=weight)) +
      geom_point(size=2, alpha=1) +
      scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(sc.weights$weight),max(sc.weights$weight))) +
      ggtitle(glue("{sc$region.disp} Cluster {sc$cluster}\nIC{sc$ic.number}")) +
      theme_few() + theme(strip.text.x=element_text(size=20), strip.text.y=element_text(size=14))
  }

  key.str <- digest(c(sc$exp.label, sc$cluster, sc$ic.number))

  renderCacheImage(tsne.plot, glue("tsne_IC_{key.str}"), 500, 500, progress=progress)
}, deleteFile=FALSE)

# display all of the cached images associated with the current clusters
output$ic.grid <- renderUI({
  # TODO: if IC selected then show a single IC on left and gene table on right, else...
  if (nrow(selected.components())>0) {
    div(style="width: 100% ; height: 500px", 
        tagList(fluidRow(column(6, imageOutput("tsne.selected.component")),
                         column(6, tabsetPanel(type="pills",
                                               tabPanel("Top Genes", DT::dataTableOutput("dt.IC.genes")),
                                               tabPanel("Subclusters", DT::dataTableOutput("dt.IC.subclusters")))))))
  } else {
    # display all ICs for the filtered clusters
  
    fn <- dlply(component.cluster(), .(c.id), function(cc) {
      file.path(basename(cache.dir),"ic",list.files(file.path(cache.dir, "ic"), glue("{cc$exp.label}_{cc$cluster}_{input$opt.components}_[0-9]+_[0-9]+.png")))
    }) %>% unlist

    div(class="img-center", style="width: 100% ; height: 500px", tagList(lapply(fn, function(f) img(src=f))))
  }
})

output$dt.components <- DT::renderDataTable( {
  req(nrow(clusters.selected.components())>0)
  components.tbl <- clusters.selected.components() %>%
    mutate(status=ifelse(status=='Real','Biological',as.character(status))) %>%
    mutate(Loadings=glue("<img height='75' width='250' src='cache/ic/ic_{exp.label}_{cluster}_IC{ic.number}_250_75.png'/>")) %>%
    dplyr::select(IC=ic.number, Region=region.disp, Cluster=cluster.disp, Class=cell_class, Status=status, "Annotation Notes"=hypothesized_common_name, Loadings, "Anatomical Notes"=anatomical_region)

  # Hack. #42 selected row is not valid on initialization. save bookmarked row and restore the first time that table is ready
  print(head(components.tbl))
  write.log('delayed.dt.components_rows_selected=',delayed.dt.components_rows_selected)
  selected.row <- delayed.dt.components_rows_selected
  if (!is.null(delayed.dt.components_rows_selected)) {
    delayed.dt.components_rows_selected <<- NULL
  }
  
  DT::datatable(components.tbl,
                rownames=FALSE, escape=FALSE,
                filter="bottom", 
                selection=list(mode="single", selected=selected.row, target='row'),
                options=list(dom="tp", pageLength=50))
})

output$dt.components.heading <- renderUI({
  if (nrow(component.cluster())>0) {
    fluidRow(
      column(10, {
        if (nrow(component.cluster()) < 10) {
          # print exp.labelA c1+c2, exp.labelB c1+c4
          clusters.by.exp <- dlply(component.cluster(), .(region.disp), function(df) df$cluster)
          tags$h4(paste("ICs for Clusters:", 
                        paste(mapply(function(exp.label, cluster) {
                          paste(exp.label, paste0(cluster, collapse=','))
                        }, names(clusters.by.exp), clusters.by.exp), collapse='; ')))
        } else {
          tags$h4(glue("ICs for {nrow(component.cluster())} clusters"))
        }
      }),
      column(2, 
             selectInput("opt.components", "Show", choices=c("Biological"='real','All'='all')))
    )
  } else {
    tags$p()
  }
})

