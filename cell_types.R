#####################################################################################################
# cell_types refers to all the different cluster and subclusters associated with different experiments
#
# There are three parts:
#  1. The dynamic UI elements for the user to choose
#  2. Filtering on the user selection and returning reactive tables, primarily clusters.selected() 
#  3. Table display


#####################################################################################################
# Cell Type Filter options used to generate the selectizeInputs
#
# the *.options reactive expressions return the corresponding options for that variable, excluding 
# rows from cell.types that are filtered out by the remaining variables.
tissue.options <- reactive({
log.reactive("fn: tissue.options")
  user.selection.changed()
  tissue.opts <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='tissue')
  return(sort(unique(tissue.opts$region.disp)))
  # tissue.vals <- setNames(tissue.opts$exp.label, tissue.opts$region.disp)
  # sort(tissue.vals[!duplicated(tissue.vals)])
})

cell.class.options <- reactive({
log.reactive("fn: cell.class.options")
  user.selection.changed()
  cell.class.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.class')
  return(sort(unique(cell.class.opt$class.disp)))
  # cell.class.vals <- setNames(cell.class.opt$class, cell.class.opt$class.disp)
  # sort(cell.class.vals[!duplicated(cell.class.vals)])
})

cell.cluster.options <- reactive({
log.reactive("fn: cell.cluster.options")
  user.selection.changed()
  cell.cluster.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.cluster')
  return(sort(unique(cell.cluster.opt$cluster.disp)))
  # cell.cluster.vals <- setNames(cell.cluster.opt$cluster, cell.cluster.opt$cluster.disp)
  # sort(cell.cluster.vals[!duplicated(cell.cluster.vals)])
})

# A cell type is 1:1 with a cell subcluster
cell.type.options <- reactive({
log.reactive("fn: cell.type.options")
  user.selection.changed()
  cell.type.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.type')
  return(sort(unique(cell.type.opt$subcluster.disp)))
  # cell.type.vals <- setNames(cell.type.opt$subcluster, cell.type.opt$subcluster.disp)
  # cell.type.vals[!duplicated(cell.type.vals)] %>% sort
})

#####################################################################################################
# OUTPUT - user filter selections on cell types
output$region <- renderUI({
  selectizeInput("tissue", "Region", choices=tissue.options(), selected=input$tissue, multiple=TRUE)
})

output$cell.class <- renderUI({
  selectizeInput("cell.class", "Cell Class", choices=cell.class.options(), selected=input$cell.class, multiple=TRUE)
})

output$cell.cluster <- renderUI({
  selectizeInput("cell.cluster", "Cell Cluster", choices=cell.cluster.options(), selected=input$cell.cluster, multiple=TRUE)
})

output$cell.type <- renderUI({
  selectizeInput("cell.type", "Cell Type", choices=cell.type.options(), selected=input$cell.type, multiple=TRUE)
})

# output$class.marker <- renderUI({
#   selectizeInput("class.marker", "Curated Class Markers", choices=class.marker.options(), selected=input$class.marker, multiple=TRUE)
# })
# 
# output$type.marker <- renderUI({
#   selectizeInput("type.marker", "Curated Type Markers", choices=type.marker.options(), selected=input$type.marker, multiple=TRUE)
# })

#####################################################################################################
#
# Filtering on cell types
#
# the name of the variable is the name in the input list
# the value is the name of the column in cell.types
#cell.types.filter.opts <- c(tissue='exp.label', cell.class='class', cell.cluster='cluster', cell.type='subcluster')
cell.types.filter.opts <- c(tissue='region.disp', cell.class='class.disp', cell.cluster='cluster.disp', cell.type='subcluster.disp')

# generate query on df using filter.opts 
# This function is called repeatedly to either return the user selected set
# or to determine what are the remaining options for a particular multichoice value
# If excl is specified, then those fields are not included.
# if only is specified, then only those fields are included in the query.
filter.df <- function(df, filter.opts, excl=NULL, only=NULL) {
  if (!is.null(excl)) {
    filter.opts <- filter.opts[!(names(filter.opts) %in% excl)]
  }
  if (!is.null(only)) {
    filter.opts <- filter.opts[names(filter.opts) %in% only]
  }
  
  # create a list from input with only the named filters and augment the 1:N values for markers
  if (class(input)!='list') {
    user.input <- isolate(reactiveValuesToList(input))
  } else {
    user.input <- input
  }
  
  # replace the value of class.marker with the corresponding multi-value fields 
  # if (!is.null(user.input$class.marker)) {
  #   user.input$class.marker <- unique(filter(class.markers, class_marker %in% user.input$class.marker))$class_markers
  # }
  # if (!is.null(user.input$type.marker)) {
  #   user.input$type.marker <- unique(filter(type.markers, type_marker %in% user.input$type.marker))$type_markers
  # }
  
  # there are different ways of doing this, but this is the easiest, particularly because
  # I'm piecing together an unknown number of different filter_ calls. 
  # See https://stackoverflow.com/a/27197858/86228 among other refs. 
  # The "problem" is that this hides the reactive input from triggering on change
  expr <- paste0('(',paste0('is.null(user.input$',names(filter.opts), ') | ',filter.opts, ' %in% user.input$', names(filter.opts), collapse=') & ('),')')
  filter_(df, expr)
}

# this dumby reactive is needed because filter.df hides the reactive variables  
user.selection.changed <- reactive({
log.reactive("fn: user.selection.changed")
  # make reactive on these user inputs
  input$tissue; input$cell.class; input$cell.cluster; input$cell.type; input$class.marker; input$type.marker
  input$use.common.name; input$opt.cluster.disp
  # FIXME: this always returns TRUE, which is fine because it triggers dependent reactives, but
  # it would be nice to be able to write if (!user.selection.changed()).
  TRUE
})

# add display labels as per user options
subclusters.labeled <- reactive({
log.reactive("fn: subclusters.labeled")
  inner_join(cell.types, region.names(), by='exp.label') %>% 
    inner_join(cluster.names(), by=c('exp.label','cluster')) %>%
    inner_join(subcluster.names(), by=c('exp.label','subcluster')) %>%
    mutate(class.disp=class)
})

limit.by.components <- function(df, comps) {
  if (nrow(comps)>0) {
    filter(df, cluster==first(comps$cluster) & exp.label==first(comps$exp.label))
  } else {
    df
  }
}

# adds columns for each user gene containing 
gene.cols <- function(df, kind) {
  lapply(input$user.genes, function(g) {
    g.fn <- glue("{prep.dir}/markers/genes/{kind}/{g}.diffexp.RDS")
    if (file.exists(g.fn)) {
      g.diffexp <- readRDS(g.fn) %>% select(exp.label, cx, log.target.u, pval)

      names(g.diffexp)[3] <- paste0(g,'.',names(g.diffexp[3]))
      names(g.diffexp)[4] <- paste0(g,'.',names(g.diffexp[4]))
      
      by.names <- c('exp.label','cx') %>% setNames(c('exp.label',kind))
      df <<- left_join(df, g.diffexp, by=by.names)  
    }
  })
  df
}

# returns a tibble of only the rows that match the user selection
subclusters.selected_ <- reactive({
  log.reactive("fn: subclusters.selected_")
  user.selection.changed()
  
  filter.df(subclusters.labeled(), cell.types.filter.opts)
})

subclusters.selected__ <- reactive({
  limit.by.components(subclusters.selected_(), selected.components())
})


subclusters.selected <- reactive({
  gene.cols(subclusters.selected__(),'subcluster')
})

# returns a smaller tibble with a row per cluster from subclusters.selected()
clusters.selected <- reactive({
log.reactive("fn: clusters.selected")
  user.selection.changed()
  select(subclusters.selected__(), region.disp, class.disp, cluster.disp, exp.label, cluster) %>% unique %>%
    gene.cols('cluster')
})

regions.selected <- reactive({
log.reactive("fn: regions.selected")
  user.selection.changed()
  select(subclusters.selected(), region.disp, exp.label) %>% unique 
})

# # The markers are indirect values
# class.marker.options <- reactive({
#   user.selection.changed()
#   (filter.df(cell.types, cell.types.filter.opts, excl='class.marker') %>%
#      inner_join(class.markers, by='class'))$class_marker.y %>% na.omit() %>% sort
# })
# 
# type.marker.options <- reactive({
#   user.selection.changed()
#   (filter.df(cell.types, cell.types.filter.opts, excl='type.marker') %>%
#      inner_join(type.markers, by='full_name'))$type_marker.y %>% na.omit() %>% sort
# })



#####################################################################################################
# OUTPUT - Cluster/Subcluster Tables
output$dt.clusters <- DT::renderDataTable({
  ct <- clusters.selected()
  col.idx <- c(
    which(names(ct) %in% c('region.disp','class.disp','cluster.disp')),
    lapply(input$user.genes, function(g) grep(paste0('^',g,'\\.'), names(ct))) %>% unlist
  )
  ct <- ct[,col.idx]
  if (ncol(ct)>3) ct <- ct[rev(order(ct[[4]])),]  # sort by first gene's transcript amount
  
  colnames <- c('Region','Class','Cluster',
                 lapply(input$user.genes, function(g) paste(g,c('Amount','P-Val'))) %>% unlist)
  
  DT::datatable(ct, 
                rownames=FALSE,
                selection="single",
                colnames=colnames,
                options=list(dom='t', paging=FALSE)    )
})

observeEvent(input$dt.clusters_rows_selected, {
  updateSelectInput(session, 'current.cluster', selected=input$dt.clusters_rows_selected)
})

output$dt.subclusters <- DT::renderDataTable({
  
  ct <- subclusters.selected()
  col.idx <- c(
    which(names(ct) %in% c('region.disp','class.disp','cluster.disp','subcluster.disp')),
    lapply(input$user.genes, function(g) grep(paste0('^',g,'\\.'), names(ct))) %>% unlist
  )
  ct <- ct[,col.idx]
  if (ncol(ct)>4) ct <- ct[rev(order(ct[[5]])),]
  
  colnames <- c('Region','Class','Cluster','Sub-Cluster',
                lapply(input$user.genes, function(g) paste(g,c('Amount','P-Val'))) %>% unlist)
  DT::datatable(ct, 
                rownames=FALSE,
                selection="single",
                colnames=colnames,
                options=list(dom='t', paging=FALSE)    )
})

observeEvent(input$dt.subclusters_rows_selected, {
  updateSelectInput(session, 'current.subcluster', selected=input$dt.subclusters_rows_selected)
})


