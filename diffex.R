# Displaying differentially expressed genes


# for the selected marker (from the differentially expressed gene list)
# return the cells, expression levels, global X, Y coordinates in the comparison region
expr.xy <- reactive({
log.reactive("fn: expr.xy")
  ddply(cluster.markers.selected(), .(GENE), function(g) {
    ddply(regions.selected(), .(exp.label), function(r) {
      expr.fn <- glue("expr/{first(r$exp.label)}/gene/{first(g$GENE)}.RDS")
      if (file.exists(expr.fn)) {
        readRDS(expr.fn) %>% mutate(exp.label=r$exp.label) %>% inner_join(region.names(), by='exp.label')
      } else {
        tibble()
      }
    })
  })
})

# for the selected subcluster markers 
# return the cells, expression levels, global X, Y limited to selected subclusters
expr.subcluster.xy <- reactive({
log.reactive("fn: expr.subcluster.xy")
  ddply(subcluster.markers.selected(), .(GENE), function(g) {
    ddply(regions.selected(), .(exp.label), function(r) {
      expr.fn <- glue("expr/{first(r$exp.label)}/gene/{first(g$GENE)}.RDS")
      if (file.exists(expr.fn)) {
        readRDS(expr.fn) %>% mutate(exp.label=r$exp.label) %>% inner_join(region.names(), by='exp.label')
      } else {
        tibble()
      }
    })
  })
})

# for the selected subcluster markers
# return the cells, expression levels and local X, Y limited to selected subclusters
expr.subcluster.local.xy <- reactive({
log.reactive("fn: expr.subcluster.local.xy")
  if (nrow(clusters.selected()) > MAX_REGIONS) {
    log("Error: Too many facets to display. Not including individual cells")
    tibble()
  } else {
    xy <- expr.subcluster.xy()
    if (nrow(xy)>0) {
      inner_join(xy %>% select(-V1, -V2, -region.disp),            # replace x,y with local ones
                 local.xy.selected(), by=c('cell','exp.label')) %>% as_tibble()
    } else {
      xy
    }
  }
})
