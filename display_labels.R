# reactive utility functions referencing global experiment data

# returns a tibble containing the currently configured display form of the region/experiment
# [ exp.label, region.disp ]
region.names <- reactive({
log.reactive("fn: region.names")
  if (input$opt.region.disp=='region') {
    tibble(exp.label=experiments$exp.label, region.disp=experiments$exp.title, region.abbrev=experiments$exp.abbrev)
  } else {
    tibble(exp.label=experiments$exp.label, region.disp=experiments$exp.label, region.abbrev=experiments$exp.abbrev)
  }
}) 

# returns a tibble containing the currently configured display form of the cluster
# [ exp.label, cluster, cluster.disp ]
cluster.names <- reactive({
log.reactive("fn: cluster.names")
  if (input$opt.cluster.disp=='numbers') {
    tibble(exp.label=cluster.names_$exp.label, cluster=cluster.names_$cluster, cluster.disp=cluster.names_$cluster)
  } else { # annotated or all
    df <- tibble(exp.label=cluster.names_$exp.label, cluster=cluster.names_$cluster, cluster.disp=cluster.names_$cluster_name)
    if (input$opt.cluster.disp=='all') {
      mutate(df, cluster.disp=sprintf("%s [#%s]", cluster.disp, cluster))
    } else
      df
  }
})


# returns a tibble for labeling clusters in plots. This is similar, but slightly different than cluster.names().
# The disp result is whatever cluster.names() returns as cluster.disp (which could be a number).
# The number result is the cluster.
# And none returns NA to inhibit plotting
cluster.labels <- reactive({
  if (input$opt.plot.label=='disp') {
    cluster.names()
  } else if (input$opt.plot.label=='number') {
    mutate(cluster.names(), cluster.disp=cluster)
  } else if (input$opt.plot.label=='none') {
    mutate(cluster.names(), cluster.disp=NA)
  } else { stop("Unknown opt.plot.label") }
})

# returns a tibble containing the currently configured display form of the subcluster
# [ exp.label, subcluster, subcluster.disp ]
subcluster.names <- reactive({
log.reactive("fn: subcluster.names")
  if (input$opt.cluster.disp=='numbers') {
    tibble(exp.label=subcluster.names_$exp.label, subcluster=subcluster.names_$subcluster, subcluster.disp=subcluster.names_$subcluster)
  } else { # annotated or all

    if (input$use.common.name)
      df <- tibble(exp.label=subcluster.names_$exp.label, subcluster=subcluster.names_$subcluster, subcluster.disp=subcluster.names_$subcluster_name)
    else
      df <- tibble(exp.label=subcluster.names_$exp.label, subcluster=subcluster.names_$subcluster, subcluster.disp=subcluster.names_$full_name)
    
    if (input$opt.cluster.disp=='all') {
      mutate(df, subcluster.disp=sprintf("%s [#%s]", subcluster.disp, subcluster))
    } else
      df
  }
})

subcluster.labels <- reactive({
  if (input$opt.plot.label=='disp') {
    subcluster.names()
  } else if (input$opt.plot.label=='number') {
    mutate(subcluster.names(), subcluster.disp=subcluster)
  } else if (input$opt.plot.label=='none') {
    mutate(subcluster.names(), subcluster.disp=NA)
  } else { stop("Unknown opt.plot.label") }
})

# to allow lookups of cluster from subcluster or vice-versa
all.subclusters <- reactive({
log.reactive("fn: all.subclusters")
  with(cell.types, tibble(exp.label=exp.label, cluster=cluster, subcluster=subcluster))
})

# a tibble with cluster and subcluster combined
# [ exp.label, cx, cx.disp ]
cx.names <- reactive({
log.reactive("fn: cx.names")
  rbind(dplyr::rename(cluster.names(), cx=cluster, cx.disp=cluster.disp), dplyr::rename(subcluster.names(), cx=subcluster, cx.disp=subcluster.disp)) %>%
    mutate(cx=as.character(cx))
})

