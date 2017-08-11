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
# The user can choose Fold Change, Percent Present or BOTH or EITHER
# This disables/enables fold change and percent present depending on the user option
observeEvent(input$expr.filter.opt, {
  if (input$expr.filter.opt=='pp') {
    shinyjs::disable("fold.change")
  } else {
    shinyjs::enable("fold.change")
  }
  
  pct.present.widgets <- c('min.pct.within','max.pct.without')
  if (input$expr.filter.opt=='fc') {
    lapply(pct.present.widgets, function(w) shinyjs::disable(w) )
  } else {
    lapply(pct.present.widgets, function(w) shinyjs::enable(w) )
  }      
})

#####################################################################################################
# Gene Expression Filter

# mark those genes that meet expression criteria
diff.exp.filter <- function(df) {
  df$fold.change.thresh <- df$FOLDch >= input$fold.change
  df$pct.thresh <- df$pct.1 >= input$min.pct.within & df$pct.2 <= input$max.pct.without
  
  # either, both, fc, pp
  if (input$expr.filter.opt == 'either')
    mutate(df, expr.pass = fold.change.thresh | pct.thresh)
  else if (input$expr.filter.opt == 'both')
    mutate(df, expr.pass = fold.change.thresh & pct.thresh)
  else if (input$expr.filter.opt == 'fc')
    mutate(df, expr.pass = fold.change.thresh)
  else if (input$expr.filter.opt == 'pp')
    mutate(df, expr.pass = pct.thresh)
  else
    stop(glue("Wrong expr.filter.opt {input$expr.filter.opt}"))
}


# returns the filename to read the marker data, depending on what kind of comparison the user wants
# kind: cluster or subcluster
# comp.type: global, pairwise
# region: exp.label
markers.fn <- function(kind, comp.type, region, cluster=NULL) {
  if (comp.type=='global') {
    fn <- glue("{prep.dir}/markers/{kind}Markers.RDS")
  } else {
    fn <- glue("{prep.dir}/markers/{region}/{cluster}.{kind}.pairwise.markers.RDS")
  }
  fn
}

# returns all comparisons or pairwise for the current target and comparison
# target: a row from clusters.selected() or subclusters.selected()
# comparison: a comparison row from clusters.selected() or subclusters.selected() or the word "global"
# kind: "cluster" or "subcluster".
#
# If "global", then returns the set of differentially expressed genes for the region of the target
all.markers <- function(target, comparison, kind) {
  if (comparison[[kind]]=='global') {
    filter(readRDS(markers.fn(kind,'global',target$exp.label)), cluster==target[[kind]] & exp.label==target$exp.label)
  } else {
    filter(readRDS(markers.fn(kind,'pairwise',target$exp.label, target[[kind]])), cluster==target[[kind]] & other.cluster==comparison[[kind]])
  }
}

filtered.markers <- function(target, comparison, kind) {
  df <- all.markers(target, comparison, kind)
  
  df$user.selected <- df$GENE %in% input$user.genes
  
  # adds a column, expr.pass if it meets expression criteria
  df <- diff.exp.filter(df)
  
  # filter those that meet exp criteria or are user selected. remove non-matching tissues
  df <- filter(df, user.selected | expr.pass)

  # add display columns
  df <- inner_join(df, region.names(), by=c('exp.label'))

  {
    if (kind=='cluster') {
      if (comparison[[kind]]=='global')
        inner_join(df, cluster.names(), by=c('exp.label','cluster'))
      else
        inner_join(df, cluster.names(), by=c('exp.label','cluster')) %>%
        inner_join(dplyr::rename(cluster.names(), other.cluster='cluster', other.cluster.disp='cluster.disp'),
                   by=c('exp.label','other.cluster'))
    } else { # subcluster, but the columns are still labeled cluster and other.cluster
      if (comparison[[kind]]=='global')
        inner_join(df, subcluster.names(), by=c('exp.label',cluster='subcluster'))
      else
        inner_join(df, subcluster.names(), by=c('exp.label',cluster='subcluster')) %>%
        inner_join(dplyr::rename(subcluster.names(), other.cluster='subcluster', other.cluster.disp='subcluster.disp'),
                   by=c('exp.label','other.cluster'))
    }
  } %>% arrange(-user.selected,-expr.pass, -pct.1)
}

# Returns a table of the markers associated with the filtered clusters, 
# filtered according to expression criteria and user genes
cluster.markers <- reactive({
log.reactive("fn: cluster.markers")
  req(current.cluster.i())
  filtered.markers(current.cluster(), comparison.cluster(), 'cluster')
})


subcluster.markers <- reactive({
log.reactive("fn: subcluster.markers")
  req(current.subcluster.i())
  filtered.markers(current.subcluster(), comparison.subcluster(), 'subcluster')
})

cluster.markers.selected <- reactive({
log.reactive("fn: cluster.markers.selected")
  if (is.null(input$dt.cluster.markers_rows_selected)) {
    tibble(GENE=character(0))
  } else {
    cluster.markers()[input$dt.cluster.markers_rows_selected,]
  }
})

subcluster.markers.selected <- reactive({
log.reactive("fn: subcluster.markers.selected")
  if (is.null(input$dt.subcluster.markers_rows_selected)) {
    tibble(GENE=character(0))
  } else {
    subcluster.markers()[input$dt.subcluster.markers_rows_selected,]
  }
})

#####################################################################################################
# OUTPUT - Markers tables

dt.markers <- function(mrkrs, other) {
  if (other=='global')
    mrkrs$other.cluster.disp <- mrkrs$exp.label
  
  DT::datatable(select(mrkrs, GENE, FOLDch, pct.1, pct.2, row.highlight),
                rownames = FALSE,
                selection="multiple",
                colnames = c('Gene','Fold Change','Target Present','Other Present','row.highlight'),
                options=list(dom="tp", pageLength=50, columnDefs = list(list(visible=FALSE, targets=4)))) %>% 
    DT::formatStyle('row.highlight', target='row', 
                    backgroundColor = DT::styleEqual(c(0,1,2), c('pink','lightgreen','white'))) %>%
    DT::formatSignif('FOLDch', 3) %>% DT::formatPercentage('pct.1', 1) %>% DT::formatPercentage('pct.2', 1)
}

output$dt.cluster.markers <- DT::renderDataTable( {
  mrkrs <- mutate(cluster.markers(), row.highlight=ifelse(user.selected,ifelse(expr.pass,1,0),2)) 
  dt.markers(mrkrs, comparison.cluster()$cluster)
})

output$dt.cluster.markers.dl <- downloadHandler(filename="cluster-markers.csv", 
                                                content= function(file) {
                                                  write.csv(cluster.markers(), file=file)
                                                })


output$dt.subcluster.markers <- DT::renderDataTable( {
  mrkrs <- mutate(subcluster.markers(), row.highlight=ifelse(user.selected,ifelse(expr.pass,1,0),2))
  dt.markers(mrkrs, comparison.subcluster()$subcluster)
})

output$dt.subcluster.markers.dl <- downloadHandler(filename="subcluster-markers.csv", 
                                                   content= function(file) {
                                                     write.csv(subcluster.markers(), file=file)
                                                   })
