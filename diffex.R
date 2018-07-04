# Displaying differentially expressed genes


# for the selected marker (from the differentially expressed gene list) or a user-supplied gene(s),
# return the cells, expression levels, global X, Y coordinates in the comparison region
expr.xy <- reactive({
  ldply(user.genes(), function(g) {
    ddply(regions.selected(), .(exp.label), function(r) {
      expr.fn <- glue("{prep.dir}/expr/{first(r$exp.label)}/gene/{g}.RDS")
      if (file.exists(expr.fn)) {
        readRDS(expr.fn) %>% mutate(exp.label=r$exp.label) %>% inner_join(region.names(), by='exp.label')
      } else {
        tibble()
      }
    })
  })
})

# for the selected subcluster markers or user-supplied gene(s),
# return the cells, expression levels, global X, Y limited to selected subclusters
## expr.subcluster.xy <- reactive({
##   stop("expr.xy and expr.subcluster.xy are now identical. merge?")
  
##   ldply(user.genes(), function(g) {
##     ddply(regions.selected(), .(exp.label), function(r) {
##       expr.fn <- glue("{prep.dir}/expr/{first(r$exp.label)}/gene/{first(g$gene)}.RDS")
##       if (file.exists(expr.fn)) {
##         readRDS(expr.fn) %>% mutate(exp.label=r$exp.label) %>% inner_join(region.names(), by='exp.label')
##       } else {
##         tibble()
##       }
##     })
##   })
## })
expr.subcluster.xy <- reactive({ expr.xy() })

# for the selected subcluster markers or user-supplied gene(s)
# return the cells, expression levels and local X, Y limited to selected subclusters
expr.subcluster.local.xy <- reactive({
  
  xy <- expr.subcluster.xy()
  if (nrow(xy)>0) {
    inner_join(xy %>% select(-V1, -V2, -region.disp),            # replace x,y with local ones
               local.xy.selected(), by=c('cell','exp.label')) %>% as_tibble()
  } else {
    xy
  }
})
