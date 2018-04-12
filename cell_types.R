#####################################################################################################
# cell_types refers to all the different cluster and subclusters associated with different experiments
#
# There are three parts:
#  1. The dynamic UI elements for the user to choose
#  2. Filtering on the user selection and returning reactive tables, primarily clusters.selected() 
#  3. Table display

filter.vals <- reactiveValues()

filter.vals.init <- function() {
  filter.vals$user.genes <- NULL
  filter.vals$tissue <- NULL
  filter.vals$cell.class <- NULL
  filter.vals$cell.cluster <- NULL
  filter.vals$cell.type <- NULL
}
filter.vals.init()

# if the user pressed the go button ("Update") or resets parameters
# then disable go button and remove dirty indictor
observeEvent(input$go, {
  filter.vals$user.genes <- isolate(input$user.genes)
  filter.vals$tissue <- isolate(input$tissue)
  filter.vals$cell.class <- isolate(input$cell.class)
  filter.vals$cell.cluster <- isolate(input$cell.cluster)
  filter.vals$cell.type <- isolate(input$cell.type)
}, ignoreInit=TRUE)

# returns true if a==b.
# blank and null are considered the same
input.cmp <- function(a, b) {
  if (is.null(a)) a <- ''
  if (is.null(b)) b <- ''
  return(all(a==b))
}

# if the user changes any of these parameters, then
# set the dirty indicator until presses go.
observe({
  if (input.cmp(filter.vals$user.genes, input$user.genes) &&
        input.cmp(filter.vals$tissue, input$tissue) &&
        input.cmp(filter.vals$cell.class, input$cell.class) &&
        input.cmp(filter.vals$cell.cluster, input$cell.cluster) &&
        input.cmp(filter.vals$cell.type, input$cell.type)) {
    write.log("Equal")
    removeCssClass("filter-params", "dirty-controls")
    disable("go")
  } else {
    write.log("Mismatch")
    write.log(glue("user.genes {ifelse(is.null(filter.vals$user.genes),'NULL',filter.vals$user.genes)} <-> {ifelse(is.null(input$user.genes),'NULL',input$user.genes)} {input.cmp(filter.vals$user.genes, input$user.genes)}"))
    write.log(glue("tissue {ifelse(is.null(filter.vals$tissue),'NULL',filter.vals$tissue)} <-> {ifelse(is.null(input$tissue),'NULL',input$tissue)} {input.cmp(filter.vals$tissue, input$tissue)}"))
    write.log(glue("cell.class {ifelse(is.null(filter.vals$cell.class),'NULL',filter.vals$cell.class)} <-> {ifelse(is.null(input$cell.class),'NULL',input$cell.class)} {input.cmp(filter.vals$cell.class, input$cell.class)}"))
    write.log(glue("cell.cluster {ifelse(is.null(filter.vals$cell.cluster),'NULL',filter.vals$cell.cluster)} <-> {ifelse(is.null(input$cell.cluster),'NULL',input$cell.cluster)} {input.cmp(filter.vals$cell.cluster, input$cell.cluster)}"))
    write.log(glue("cell.type {ifelse(is.null(filter.vals$cell.type),'NULL',filter.vals$cell.type)} <-> {ifelse(is.null(input$cell.type),'NULL',input$cell.type)} {input.cmp(filter.vals$cell.type, input$cell.type)}"))
    addCssClass("filter-params", "dirty-controls")
    enable("go")
  }
})

# set user.genes to stored value, if present
# because the options for user.genes are delay-loaded, then the initial restore
# of input will fail because the values for user.genes are missing from the choices.
# But later, with input$user.genes blank, this event will fire and javascript
# will be called (once) to reset the user.genes by forcing values.
observeEvent(input$user.genes, {
  if (!isTruthy(input$user.genes) && isTruthy(filter.vals$user.genes)) {
    js$setgenes(items=as.list(filter.vals$user.genes))
  }
    
}, priority=-100, ignoreNULL=FALSE, once=TRUE)
  
onRestore(function(state) {
  filter.vals$user.genes <- state$input$user.genes
  filter.vals$tissue <- state$input$tissue
  filter.vals$cell.class <- state$input$cell.class
  filter.vals$cell.cluster <- state$input$cell.cluster
  filter.vals$cell.type <- state$input$cell.type
})



#####################################################################################################
# Cell Type Filter options used to generate the selectizeInputs
#
# the *.options reactive expressions return the corresponding options for that variable, excluding 
# rows from cell.types that are filtered out by the remaining variables.
tissue.options <- reactive({
log.reactive("fn: tissue.options")
  user.selection.changed()
  tissue.opts <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='tissue', use.input=TRUE)
  return(sort(unique(tissue.opts$region.disp)))
  # tissue.vals <- setNames(tissue.opts$exp.label, tissue.opts$region.disp)
  # sort(tissue.vals[!duplicated(tissue.vals)])
})

cell.class.options <- reactive({
log.reactive("fn: cell.class.options")
  user.selection.changed()
  cell.class.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.class', use.input=TRUE)
  return(sort(unique(cell.class.opt$class.disp)))
  # cell.class.vals <- setNames(cell.class.opt$class, cell.class.opt$class.disp)
  # sort(cell.class.vals[!duplicated(cell.class.vals)])
})

cell.cluster.options <- reactive({
log.reactive("fn: cell.cluster.options")
  user.selection.changed()
  cell.cluster.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.cluster', use.input=TRUE)
  
  # prefix region abbrev if more than one region to distinguish, e.g. "Neuron [#2]" => "STR Neuron [#2]"
  if (length(filter.vals$tissue) != 1) {
    cco <- arrange(unique(select(cell.cluster.opt, region.abbrev, cluster.disp)), cluster.disp)
    return(setNames(cco$cluster.disp, glue("{cco$region.abbrev} {cco$cluster.disp}")))
  } else {
    return(sort(unique(cell.cluster.opt$cluster.disp)))
  }
  # cell.cluster.vals <- setNames(cell.cluster.opt$cluster, cell.cluster.opt$cluster.disp)
  # sort(cell.cluster.vals[!duplicated(cell.cluster.vals)])
})

# A cell type is 1:1 with a cell subcluster
cell.type.options <- reactive({
log.reactive("fn: cell.type.options")
  user.selection.changed()
  cell.type.opt <- filter.df(subclusters.labeled(), cell.types.filter.opts, excl='cell.type', use.input=TRUE)
  return(sort(unique(cell.type.opt$subcluster.disp)))
  # cell.type.vals <- setNames(cell.type.opt$subcluster, cell.type.opt$subcluster.disp)
  # cell.type.vals[!duplicated(cell.type.vals)] %>% sort
})

#####################################################################################################
# OUTPUT - user filter selections on cell types
output$region <- renderUI({
  selectizeInput("tissue", "Limit By Region", choices=c("Tissue / Region"="",tissue.options()), selected=input$tissue, multiple=TRUE)
})

output$cell.class <- renderUI({
  selectizeInput("cell.class", "Limit By Class", choices=c("Cell Class"="",cell.class.options()), selected=input$cell.class, multiple=TRUE)
})

output$cell.cluster <- renderUI({
  selectizeInput("cell.cluster", "Limit By Cluster", choices=c("Cell Cluster"="",cell.cluster.options()), selected=input$cell.cluster, multiple=TRUE)
})

output$cell.type <- renderUI({
  selectizeInput("cell.type", "Limit By Cell Type", choices=c("Cell Type / Subcluster"="",cell.type.options()), selected=input$cell.type, multiple=TRUE)
})


#####################################################################################################
#
# Filtering on cell types
#
# the name of the variable is the name in the input list
# the value is the name of the column in cell.types
cell.types.filter.opts <- c(tissue='region.disp', cell.class='class.disp', cell.cluster='cluster.disp', cell.type='subcluster.disp')

# generate query on df using filter.opts 
# This function is called repeatedly to either return the user selected set
# or to determine what are the remaining options for a particular multichoice value
# If excl is specified, then those fields are not included.
# if only is specified, then only those fields are included in the query.
filter.df <- function(df, filter.opts, excl=NULL, only=NULL, use.input=getOption("dropviz.use.input", default=FALSE)) {
  if (!is.null(excl)) {
    filter.opts <- filter.opts[!(names(filter.opts) %in% excl)]
  }
  if (!is.null(only)) {
    filter.opts <- filter.opts[names(filter.opts) %in% only]
  }
  
  # create a list from input or filter.vals
  if (use.input) {
    if (class(input)!='list') {
      user.input <- isolate(reactiveValuesToList(input))
    } else {
      user.input <- input
    }
  } else {
    if (class(filter.vals)!='list') {
      user.input <- isolate(reactiveValuesToList(filter.vals))
    } else {
      user.input <- filter.vals
    }
  }

  # lapply(names(filter.opts), function(nm) {
  #   write.log("user.input$",nm,"=",user.input[[nm]]," NULL:",is.null(user.input[[nm]]))
  # })
  
  # there are different ways of doing this, but this is the easiest, particularly because
  # I'm piecing together an unknown number of different filter_ calls. 
  # See https://stackoverflow.com/a/27197858/86228 among other refs. 
  # The "problem" is that this hides the reactive input from triggering on change
  expr <- paste0('(',paste0('is.null(user.input$',names(filter.opts), ') | ',filter.opts, ' %in% user.input$', names(filter.opts), collapse=') & ('),')')
  # write.log("expr: ",expr)
  filter_(df, expr)
}

# this dumby reactive is needed because filter.df hides the reactive variables  
user.selection.changed <- reactive({
  # make reactive on these user inputs
  input$tissue; input$cell.class; input$cell.cluster; input$cell.type; input$class.marker; input$type.marker
  input$use.common.name; input$opt.cluster.disp
  # FIXME: this always returns TRUE, which is fine because it triggers dependent reactives, but
  # it would be nice to be able to write if (!user.selection.changed()).
  TRUE
})

filter.vals.changed <- reactive({
  filter.vals$tissue; filter.vals$cell.class; filter.vals$cell.cluster; filter.vals$cell.type; input$class.marker; input$type.marker
  input$use.common.name; input$opt.cluster.disp
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
  # FIXME: further filtering based on component selection is obsolete
  if (FALSE && nrow(comps)>0) {
    filter(df, cluster==first(comps$cluster) & exp.label==first(comps$exp.label))
  } else {
    df
  }
}

# returns a list of gene symbols from filter.vals$user.genes
user.manual.genes <- reactive({
  if (isTruthy(filter.vals$user.genes)) {
    genes <- gene.symbol(filter.vals$user.genes)
    if (any(genes=='')) {
      missing <- filter.vals$user.genes[genes=='']
      genes <- genes[genes!='']
      missing.txt <- glue("Ignoring unknown gene symbols {paste(missing,collapse=', ')}")
      if (!is.null(getDefaultReactiveDomain()))
        showNotification(missing.txt, duration=15, type='warning')
      write.log(missing.txt)
    }
    genes
  } else {
    filter.vals$user.genes
  }
})

# return a list of genes combining user.manual.genes and any rows selected from markers tables
user.genes <- reactive({
  ug <- user.manual.genes()
  if (isTruthy(current.cluster.i())) {
    ug <- c(ug, cluster.markers.selected()$gene)
  }
  if (isTruthy(current.subcluster.i())) {
    ug <- c(ug, subcluster.markers.selected()$gene) 
  }
  unique(ug)
})

# adds columns for each user gene and sorts on the first
gene.cols <- function(df, kind) {
  lapply(user.genes(), function(g) {
    g.fn <- glue("{prep.dir}/markers/genes/{kind}/{g}.diffexp.RDS")
    if (file.exists(g.fn)) {
      g.diffexp <- readRDS(g.fn) %>% select(exp.label, cx, log.target.u, pval, target.sum.per.100k, target.sum.L.per.100k, target.sum.R.per.100k)

      # prepend gene name on columns
      lapply(3:7, function(idx) names(g.diffexp)[idx] <<- paste0(g,'_',names(g.diffexp[idx])))
      
      by.names <- c('exp.label','cx') %>% setNames(c('exp.label',kind))
      df <<- left_join(df, g.diffexp, by=by.names)  
    }
  })

  if (isTruthy(user.genes()))
    df <- df[order(df[[first(grep('.log.target.u',colnames(df)))]], decreasing = TRUE),]  # sort by first gene's transcript amount
  
  df
}

# returns a tibble of only the rows that match the user selection
# HACK: FIX ME. This is a mess because of limit.by.components. Get rid of limit.by.components and remove subclusters.selected_
# This is even more of a mess because some calls to (sub)clusters.selected() cannot include the gene.cols. It probably needs
# to stay even after the components kruft is removed.
subclusters.selected_ <- reactive({
  log.reactive("fn: subclusters.selected_")
  filter.vals.changed()
  
  filter.df(subclusters.labeled(), cell.types.filter.opts)
})

subclusters.selected__ <- reactive({
  limit.by.components(subclusters.selected_(), selected.components())
})


# add gene.cols
subclusters.selected <- reactive({
  gene.cols(subclusters.selected__(),'subcluster')
})

# returns a smaller tibble with a row per cluster from subclusters.selected()
clusters.selected_ <- reactive({
  log.reactive("fn: clusters.selected_")
  filter.vals.changed()
  select(subclusters.selected__(), c.id, region.disp, region.abbrev, class.disp, cluster.disp, exp.label, cluster) %>% unique
})

# adds gene.cols
clusters.selected <- reactive({
  log.reactive("fn: clusters.selected")

  gene.cols(clusters.selected_(), 'cluster')
})

regions.selected <- reactive({
log.reactive("fn: regions.selected")
  filter.vals.changed()
  select(subclusters.selected(), region.disp, exp.label) %>% unique 
})


#####################################################################################################
# OUTPUT - Cluster/Subcluster Tables


# set significant digits for Amount and PVal columns, if there are any
setSig <- function(dt, ct.names, start.col, end.col, sig.digits) {
  if (end.col >= start.col) {
    DT::formatSignif(dt, ct.names[start.col:end.col], sig.digits)
  } else {
    dt
  }
}

dt.clusters.fmt <- reactive({
  ct <- clusters.selected()
  col.idx <- c(
    sapply(c('region.disp','class.disp','cluster.disp'),
           function(nm) which(names(ct) %in% nm)),
    sapply(user.genes(), function(g) 
      which(grepl(paste0('^',g,'_'), names(ct)) & !grepl('target.sum', names(ct)))
    ) %>% unlist
  )
  ct <- ct[,col.idx]
  
  setNames(ct,c('Region','Class','Cluster',
                lapply(user.genes(), function(g) paste(g,c('Amount','P-Val'))) %>% unlist))

})

output$dt.clusters <- DT::renderDataTable({
  
  ct <- dt.clusters.fmt()
  DT::datatable(ct, 
                rownames=FALSE,
                selection="none",
                colnames=names(ct),
                options=list(dom='t', paging=FALSE)) %>% setSig(names(ct), 4, ncol(ct), 3)
})

output$dt.clusters.dl <- downloadHandler(filename="clusters.csv",
                                         content=function(file) {
                                           write.csv(dt.clusters.fmt(), file=file, row.names=FALSE)
                                         })

## observeEvent(input$dt.clusters_rows_selected, {
##   updateSelectInput(session, 'current.cluster', selected=input$dt.clusters_rows_selected)
## })

dt.subclusters.fmt <- reactive({
  ct <- subclusters.selected_() %>% gene.cols('subcluster')  # HACK: fixme. The subclusters.selected routines is a mess due to filtering on selected component. Needs to be cleaned up.
  
  col.idx <- c(
    sapply(c('region.disp','class.disp','cluster.disp','subcluster.disp'),
           function(nm) which(names(ct) %in% nm)),
    sapply(user.genes(), function(g) 
      which(grepl(paste0('^',g,'_'), names(ct)) & !grepl('target.sum', names(ct)))
    ) %>% unlist 
  )
  ct <- ct[,col.idx]
  
  setNames(ct, c('Region','Class','Cluster','Sub-Cluster',
                 lapply(user.genes(), function(g) paste(g,c('Amount','P-Val'))) %>% unlist))

})

output$dt.subclusters <- DT::renderDataTable({
  
  ct <- dt.subclusters.fmt()
  DT::datatable(ct,
                rownames=FALSE,
                selection="none",
                colnames=names(ct),
                options=list(dom='t', paging=FALSE)) %>% setSig(names(ct), 4, ncol(ct), 3)
})

output$dt.subclusters.dl <- downloadHandler(filename="subclusters.csv",
                                            content=function(file) {
                                              write.csv(dt.subclusters.fmt(), file=file, row.names=FALSE)
                                         })
## observeEvent(input$dt.subclusters_rows_selected, {
##   updateSelectInput(session, 'current.subcluster', selected=input$dt.subclusters_rows_selected)
## })


