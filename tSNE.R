#############################################
# tSNE plot related functions

downsample <- reactive({
  log.reactive("fn: downsample")
  if (input$opt.downsampling.method != 'none') {
    input$downsampling
  } else {
    1e9
  }
})

#############################################
# Bag Data

# returns all the bag data for all experiments. (generated in prep-tSNE.R)
# FIXME: hack to exclude clusters with no subclusters
tsne.cluster.bag.data <- 
  lapply(readRDS(glue("{prep.dir}/tsne/global.clusters.bags.Rdata")), function(df) inner_join(df, unique(select(cell.types,exp.label,cluster)), by=c('exp.label','cluster')))

tsne.subcluster.bag.data <- 
  lapply(readRDS(glue("{prep.dir}/tsne/global.subclusters.bags.Rdata")), function(df) inner_join(df, unique(select(cell.types,exp.label,subcluster)), by=c('exp.label','subcluster')))

tsne.local.bag.data <- 
  lapply(readRDS(glue("{prep.dir}/tsne/local.subclusters.bags.Rdata")), function(df) inner_join(df, cell.types, by=c('exp.label','cluster','subcluster')))

rm.zero <- function(df) {
  filter(df, !(x==0 & y==0))
}

# returns the selected bag, loop or center with additional region.disp, cluster.disp and/or subcluster.disp
# normalizes cluster/subcluster to "cx"
# add cluster facet for local
polygons.labeled <- function(polygon, regions, clusters, kind, region.nms, cluster.nms, facet.nms=NULL) {
  if (kind=='cluster') {
    filter(polygon, exp.label %in% regions$exp.label) %>% inner_join(region.nms, by='exp.label') %>%
      left_join(select(clusters, exp.label, cluster) %>%
                  mutate(cx.selected=TRUE),
                by=c('exp.label','cluster')) %>%
      mutate(cx.gg=factor(ifelse(is.na(cx.selected),NA,as.character(cluster)),levels=levels(polygon$cluster))) %>%
      inner_join(cluster.nms, by=c('exp.label','cluster')) %>%
      dplyr::rename(cx=cluster, cx.disp=cluster.disp)
  } else if (kind=='subcluster') {
    filter(polygon, exp.label %in% regions$exp.label) %>% inner_join(region.nms, by='exp.label') %>%
      left_join(select(clusters, exp.label, subcluster) %>%
                  mutate(cx.selected=TRUE),
                by=c('exp.label','subcluster')) %>%
      mutate(cx.gg=factor(ifelse(is.na(cx.selected) | (x==0 & y==0),NA,as.character(subcluster)), levels=levels(polygon$subcluster))) %>%
      inner_join(cluster.nms, by=c('exp.label','subcluster')) %>%
      dplyr::rename(cx=subcluster, cx.disp=subcluster.disp)
  } else {
    stopifnot(kind=='local')
    major.clusters <- select(clusters, exp.label, cluster) %>% unique
    x <- filter(polygon, exp.label %in% regions$exp.label) %>% inner_join(region.nms, by='exp.label') %>%
      right_join(select(major.clusters, exp.label, cluster), by=c('exp.label','cluster')) %>%
      left_join(select(clusters, exp.label, subcluster) %>% mutate(cx.selected=TRUE), by=c('exp.label','subcluster')) %>%
      mutate(cx.gg=factor(ifelse(is.na(cx.selected) | (x==0 & y==0),NA,as.character(subcluster)), levels=levels(polygon$subcluster))) %>%
      inner_join(cluster.nms, by=c('exp.label','subcluster')) %>%
      inner_join(facet.nms, by=c('exp.label','cluster')) %>%
      dplyr::rename(cx=subcluster, cx.disp=subcluster.disp, facet.gg=cluster.disp) 
    
    if (input$opt.cluster.disp=='annotated')
      x$facet.gg=glue("{x$facet.gg} [#{x$cluster}]")
    
    x
  }
}


# filters tsne.cluster.bag.data for the bag data of selected experiments
global.selected.bag <- reactive({
log.reactive("fn: global.selected.bag")
  polygons.labeled(tsne.cluster.bag.data$bags, regions.selected(), clusters.selected(), 'cluster', region.names(), cluster.labels())
})
  
global.selected.loop <- reactive({
log.reactive("fn: global.selected.loop")
  polygons.labeled(tsne.cluster.bag.data$loops, regions.selected(), clusters.selected(), 'cluster', region.names(), cluster.labels())
})

global.selected.center <- reactive({
log.reactive("fn: global.selected.center")
  polygons.labeled(tsne.cluster.bag.data$centers, regions.selected(), clusters.selected(), 'cluster', region.names(), cluster.labels())
})

global.selected.sub.bag <- reactive({
log.reactive("fn: global.selected.sub.bag")
  polygons.labeled(tsne.subcluster.bag.data$bags, regions.selected(), subclusters.selected(), 'subcluster', region.names(), subcluster.labels())
})

global.selected.sub.loop <- reactive({
log.reactive("fn: global.selected.sub.loop")
  polygons.labeled(tsne.subcluster.bag.data$loops, regions.selected(), subclusters.selected(), 'subcluster', region.names(), subcluster.labels())
})

global.selected.sub.center <- reactive({
log.reactive("fn: global.selected.sub.center")
  polygons.labeled(tsne.subcluster.bag.data$centers, regions.selected(), subclusters.selected(), 'subcluster', region.names(), subcluster.labels())
})

local.selected.sub.bag <- reactive({
log.reactive("fn: local.selected.sub.bag")
  polygons.labeled(tsne.local.bag.data$bags, regions.selected(), subclusters.selected(), 'local', region.names(), subcluster.labels(), cluster.names())
})

local.selected.sub.loop <- reactive({
log.reactive("fn: local.selected.sub.loop")
  polygons.labeled(tsne.local.bag.data$loops, regions.selected(), subclusters.selected(), 'local', region.names(), subcluster.labels(), cluster.names())
})

local.selected.sub.center <- reactive({
log.reactive("fn: local.selected.sub.center")
  polygons.labeled(tsne.local.bag.data$centers, regions.selected(), subclusters.selected(), 'local', region.names(), subcluster.labels(), cluster.names())
})


################################################################
# XY coordinates for individual cells

# the coordinates for current regions
global.xy <- reactive({
log.reactive("fn: global.xy")
  ldply(regions.selected()$exp.label, function(region) {
    # read the global XYs 
    df.region <- 
      ddply(readRDS(glue("{prep.dir}/tsne/{region}/global.xy.RDS")), .(cluster), function(df) {
        if (input$opt.downsampling.method=='cluster') {
          sample_n(df, min(nrow(df),downsample()))
        } else {
          df
        }
      })
    (if (input$opt.downsampling.method=='uniform') {
      sample_n(df.region, min(nrow(df.region), downsample()))
    } else {
      df.region
    }) %>% mutate(exp.label=factor(region,levels=levels(region.names()$exp.label)))
  }) %>% as_tibble
})

# the global.xy limited to clusters.selected()
# may set X and Y to NA because there are no cells assigned to cluster
global.xy.cluster.selected <- reactive({
  log.reactive("fn: global.xy.cluster.selected")
  left_join(select(clusters.selected(),exp.label, cluster), global.xy(), by=c('exp.label','cluster')) %>%
    inner_join(region.names(), by='exp.label') %>%
    dplyr::rename(cx=cluster)
})

# All the xy for the the current subclusters.selected() and all non-assigned cells for the clusters.selected().
# Unassigned are displayed in grey. All cells are shown in the cluster because only a subset of the assigned cells
# were used in the subclustering.
global.xy.subcluster.selected <- reactive({
  log.reactive("fn: global.xy.subcluster.selected")
  bind_rows(left_join(select(subclusters.selected(), exp.label, subcluster), global.xy(), by=c('exp.label','subcluster')),
            (left_join(select(clusters.selected(), exp.label, cluster), global.xy(), by=c('exp.label','cluster')) %>% filter(is.na(subcluster)))) %>%
    inner_join(region.names(), by='exp.label') %>% 
    select(-cluster) %>% dplyr::rename(cx='subcluster') 
})

# local XY of all selected clusters
local.xy <- reactive ({
  
  ddply(clusters.selected(), .(exp.label), function(r.cx) {
    r.xy <-
      ddply(r.cx, .(cluster), function(df) {
        local.xy.fn <- glue("{prep.dir}/tsne/{first(df$exp.label)}/cluster{first(df$cluster)}.xy.RDS")
        if (file.exists(local.xy.fn)) {
          df <- readRDS(local.xy.fn)
          if (input$opt.downsampling.method=='cluster') {
            sample_n(df, min(nrow(df), downsample()))
          } else {
            df
          }
        } else {
          warning("Missing ",local.xy.fn)
          data.frame()
        }
      })
    if (input$opt.downsampling.method=='uniform') {
      sample_n(r.xy, min(nrow(r.xy), downsample()))
    } else {
      r.xy
    }
  }) %>% as_tibble
})

# local XY of all selected subclusters - with display names added
local.xy.selected <- reactive({
log.reactive("fn: local.xy.selected")
  x <- left_join(select(subclusters.selected(),exp.label, subcluster), local.xy(), by=c('exp.label','subcluster')) %>%
    inner_join(region.names(), by='exp.label') %>%
    inner_join(cluster.names(), by=c('exp.label','cluster')) %>%
    dplyr::rename(cx=subcluster, facet.gg=cluster.disp) 
  
  if (input$opt.cluster.disp=='annotated')
    x$facet.gg=glue("{x$facet.gg} [#{x$cluster}]")
  
  x
})

opt.tx <- reactive({ isTruthy(user.genes()) && !is.null(input$opt.tx) })
tx.cells <- reactive({ input$opt.tx.cells })
tx.alpha <- reactive({ opt.tx() && input$opt.tx=='alpha' })
tx.heat <- reactive({ opt.tx() && (input$opt.tx=='heat') })
tx.scale <- reactive({ input$opt.tx.scale })
tx.legend <- reactive({ if (input$opt.tx.legend) 'legend' else 'none' })
tx.facet2 <- reactive({ opt.tx() && (!input$opt.tx.sum || tx.cells()) && user.genes() > 1 })
facet2.count <- reactive({ (if ((length(user.genes())>0 && tx.cells()) || tx.facet2()) length(user.genes()) else 0) }) # FIXME: clean up logic

# returns the sum of all of the log normal transcript counts for all user.genes
psum.amounts <- function(cx, amounts) {
  amounts$total <- log(rowSums(as_tibble(lapply(as.list(amounts), function(a) exp(a)))))
  cbind(cx, alpha=amounts$total)
}

# peg the range to always start at 0 by adding a 0 value and then removing it
cut0 <- function(x, ...) {
  cut(c(x,0), ...)[1:length(x)]
}

HEAT.COLOR.N <- 9
cluster.transcript.amounts <- reactive({
  # If a gene search that includes cell expression for more than one gene, then return amounts
  # per gene. Otherwise, sum the levels across genes
  breaks <- if (tx.scale()=='fixed') seq(0,7) else HEAT.COLOR.N
  (
    if (tx.facet2()) {
      select(clusters.selected(), exp.label, cx=cluster, ends_with('_log.target.u')) %>%
        gather(gene, alpha, ends_with('_log.target.u')) %>%
        separate(gene, 'facet2.gg', sep='_', extra='drop')
    } else {
      cx <- select(clusters.selected(), exp.label, cx=cluster)
      amounts <- select(clusters.selected(), ends_with('_log.target.u'))
      psum.amounts(cx, amounts)
    }
  ) %>% mutate(heat=cut0(alpha, breaks, include.lowest=TRUE))
})

# TODO: factor into single function
subcluster.transcript.amounts <- reactive({
  breaks <- if (tx.scale()=='fixed') seq(0,7) else HEAT.COLOR.N

  if (tx.facet2()) {
    select(subclusters.selected(), exp.label, cx=subcluster, ends_with('_log.target.u')) %>%
      gather(gene, alpha, ends_with('_log.target.u')) %>%
      separate(gene, 'facet2.gg', sep='_', extra='drop') %>%
      mutate(heat=cut(alpha, breaks))
  } else {
    cx <- select(subclusters.selected(), exp.label, cluster=cluster, cx=subcluster)
    amounts <- select(subclusters.selected(), ends_with('_log.target.u'))
    pa <- psum.amounts(cx, amounts)
    mutate(pa, heat=cut(alpha, breaks)) %>% select(-cluster)
  }
})


################################################################
# Draw tSNE

alpha.na2zero <- function(df) mutate(df, alpha=ifelse(is.na(alpha),0,alpha))

# join with alpha only if lhs has data and then replace all NA alphas with zero
left_join_alpha_heat <- function(lhs, rhs) {
  if (nrow(lhs) > 0) {
    # add a facet2.gg (i.e. a gene name) for every (sub)cluster
    # then include the alpha and heat values. this is a bit of hack
    # so that the unselected clusters display in grey for gene searches
    full_join(lhs, unique(select(rhs, exp.label, facet2.gg)), by=c('exp.label')) %>%
      left_join(rhs, by=c('exp.label','facet2.gg','cx')) %>% alpha.na2zero()
  } else {
    mutate(lhs, alpha=double())
  }
}

# when all cells are plotted, then a point size of 0.5 looks good. But when we are downsampling,
# make the points larger.
CELL.MIN.SIZE <- 0.5
CELL.MIN.SAMPLE <- 1000
CELL.MAX.SIZE <- 2.5
CELL.CEIL.SAMPLE <- 20000
xy.cell.size <- reactive({
  if (input$opt.downsampling.method=='none') {
    CELL.MIN.SIZE
  } else {
    # anything over, say, 20k is 0.5
    # make 1000 (the minimum) be 2.5
    ds <- (CELL.CEIL.SAMPLE-min(downsample()-CELL.MIN.SAMPLE,CELL.CEIL.SAMPLE))/CELL.CEIL.SAMPLE  # range 1..0
    0.5 + 2*ds
  }
  
})


tsne.disp.opts <- reactive({
log.reactive("fn: tsne.disp.opts")
  c(input$opt.cluster.disp, input$use.common.name, downsample(), input$opt.downsampling.method, 
    input$opt.region.disp, input$opt.plot.label, input$use.bag.plot, input$opt.expr.size,
    input$opt.show.bags,input$opt.tx,input$opt.tx.min,input$opt.tx.sum, input$opt.tx.scale, input$opt.tx.legend,
    input$opt.tx.cells, input$opt.cell.display.type, input$opt.expr.size, input$opt.detection.thresh)
})

## returns a function that returns either (1) a ggplot object to draw a tsne plot or (2) a closure of the ggplot object and all dependencies
## 
tsne.label <- function(is.global=TRUE, show.subclusters=FALSE, show.cells=TRUE, show.bags=FALSE, diff.genes=tibble(), return.closure=FALSE) {
  function(progress=NULL) {
    stopifnot(is.global || show.subclusters) # can't show clusters on local tsne
    
    if (opt.tx() && input$opt.tx %in% c('heat')) {
      show.bags <- TRUE
      show.cells <- FALSE
    }
    if (!(show.bags || show.cells)) return(plot.text("No data to display."))
    
    # for both global.xy (the positions of each cell) and global.[sub]cluster.avg.xy (the center position of each [sub]cluster),
    # limit to selected [sub]clusters and add pretty region name
    write.log(glue("Building tsne.label - is.global={is.global} show.subclusters={show.subclusters} show.bags={show.bags} show.cells={show.cells} diff.genes={nrow(diff.genes)}"))
    xy.data <-
      (if (show.cells) {
        if (!is.null(progress)) progress$inc(0.2, detail="Reading XY data")
        if (show.subclusters) {
          if (is.global) {
            global.xy.subcluster.selected() %>% filter(!is.na(V1)) 
          } else {
            local.xy.selected()
          }
        } else {
          global.xy.cluster.selected()
        }
      } else {
        tibble()
      })

    if (show.cells && nrow(xy.data)==0) return(plot.text("No cell data to display."))
    
    if (!is.null(progress) && nrow(xy.data)>0) progress$inc(0.2, detail=glue("Read {nrow(xy.data)} cells"))
    
    # labels
    label.data <-
      (
        if (show.subclusters) {
          if (is.global) {
            filter(global.selected.sub.center(), !is.na(cx.gg))
          } else {
            mutate(local.selected.sub.center(), 
                   cx.disp=ifelse(is.na(cx.gg), 'No Data', as.character(cx.disp)))
          }
        } else {
          na.omit(global.selected.center())
        }
      )

    if (show.subclusters) {
      if (is.global) {
        bag.data <- global.selected.sub.bag()
        loop.data <- global.selected.sub.loop()
        center.data <- global.selected.sub.center()
      } else {
        bag.data <- local.selected.sub.bag()
        loop.data <- local.selected.sub.loop()
        center.data <- local.selected.sub.center()
      }
    } else {
      bag.data <- global.selected.bag()
      loop.data <- global.selected.loop()
      center.data <- global.selected.center()
    }
    if (!is.null(progress)) progress$inc(0.2, detail=glue("{nrow(center.data)} cx bag data"))
    
    # alpha is either fixed or set by transcript amounts.
    # if fixed, then alpha range is set by scale further below.
    # if user.genes are specified, then sum all of the amounts per cx and use that as the alpha for colors
    if (tx.alpha() || tx.heat()) {
      tx.cx <- (
        if (show.subclusters) {
          subcluster.transcript.amounts() 
        } else {
          cluster.transcript.amounts()
        }
      )

      label.data <- left_join_alpha_heat(label.data, tx.cx) 
      center.data <- left_join_alpha_heat(center.data, tx.cx) 
      bag.data <- left_join_alpha_heat(bag.data, tx.cx)
      loop.data <- left_join_alpha_heat(loop.data, tx.cx)
      xy.data <- left_join_alpha_heat(xy.data, tx.cx)
      
      bag.data <- mutate(bag.data, alpha=pmax(0,alpha-0.5))
      loop.data <- mutate(loop.data, alpha=pmax(0,alpha-1))
      xy.data <- mutate(xy.data, alpha=pmax(0,alpha-1))

      # only show labels where the cluster level is greater than thresh, or there's no data
      # if no data, then set heat and alpha to NA
      label.data$pass <- with(label.data, (as.integer(heat)/length(levels(heat)))>=(input$opt.tx.min/100))
      label.data <- filter(label.data, cx.disp=='No Data' | pass)
      label.data <- mutate(label.data, alpha=ifelse(cx.disp=='No Data', NA, alpha),
                           heat=factor(ifelse(cx.disp=='No Data', NA, as.character(heat)), levels=levels(heat))) # ugh, show me a cleaner way.

      if (!is.null(progress)) progress$inc(0.2, detail=glue("Computing alpha for {user.genes()}"))
    }
    
    # diff exp are row facets
    facet2.vals <- c()
    if (nrow(diff.genes)>0) {
      diff.genes$facet2.gg <- diff.genes$gene
      facet2.vals <- c(facet2.vals, unique(diff.genes$facet2.gg))
    }
    
    # labels to appear in top left of each plot to augment vertical facet label
    if (length(facet2.vals)>0) {
      min.x <- min(loop.data$x)
      max.y <- max(loop.data$y)
      facet.label.data <- tibble(x=min.x, y=max.y, facet2.gg=facet2.vals)
    } else {
      facet.label.data <- tibble(x=numeric(), y=numeric(), facet2.gg=character())
    }
    
    # define local vars for objects not in local scope.
    # p.func must depend on only local vars, which get contained in closure's environment 
    opt.expr.size <- input$opt.expr.size
    opt.show.cells <- show.cells
    opt.show.bags <- show.bags
    opt.global <- is.global
    opt.plot.label <- input$opt.plot.label
    opt.cell.display.type <- input$opt.cell.display.type
    diff.data <- (if (opt.cell.display.type=='detect') filter(diff.genes, transcripts > input$opt.detection.thresh) else diff.genes)
    opt.horiz.facet <- (nrow(diff.data)>0 && tx.cells()) || tx.facet2()
    opt.tx.cells <- tx.cells()
    opt.tx.alpha <- tx.alpha()
    opt.tx.heat <- tx.heat()
    opt.tx.scale <- tx.scale()
    opt.tx.legend <- tx.legend()
    opt.xy.cell.size <- xy.cell.size()

    if (opt.tx.scale=='gene') showNotification("Scaling Per Gene Not Yet Implemented", duration=15, type='warning')
    
    p.func <- function() {
      source("tSNE-plot.R", local=TRUE)
    }

    if (return.closure) {
      p.func
    } else {
      if (!is.null(progress)) progress$inc(0.2, message="ggplot", detail="Rendering")
      p.func()
    }
  }
}

# TRUE if user has narrowed filter/highlight selection
is.filtered <- function() {
  isTruthy(filter.vals$tissue) || isTruthy(filter.vals$cell.class) || isTruthy(filter.vals$cell.cluster) || isTruthy(filter.vals$cell.type)
}

# Set the image size to create square facets. The display is divided into at most 3 facets across when wrapping (e.g. 333x333 for a 1000px region).
# If there's a second facet, then a facet_grid display is used.
# If the horizontal facets are greater than 4, then the square facets become too small, so the width grows, too.
tsne.image.size <- function(facet1, facet2, display.width) {
  if (facet2 > 1) {
    # grid: 
    facet.wide = (if (facet1==1) display.width %/% 2 else if (facet1 > 4) display.width %/% 4 * facet1 else display.width)
    each.width <- facet.wide %/% facet1
    facet.high <- each.width * facet2
  } else {
    # wrap
    facet.wide.count <- min(3, facet1)
    facet.wide <- ifelse(facet.wide.count==1, display.width %/% 2, display.width)
    facet.high <- (facet.wide %/% facet.wide.count) * (((facet1-1) %/% facet.wide.count)+1)
  }

  write.log(glue("facet1={facet1} facet2={facet2} facet.wide={facet.wide} facet.high={facet.high}"))
  return(list(width=facet.wide, height=facet.high))
}

#########################################################
# show global tSNE with just the filtered clusters 
output$tsne.global.cluster.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())
  
  tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=FALSE, show.cells=downsample()>0, show.bags = input$opt.show.bags, diff.genes=expr.xy())
  
  region.count <- nrow(regions.selected())
  gene.count <- facet2.count()

  display.width <- img.size.round(session$clientData[[glue("output_tsne.global.cluster.label_width")]])
  img.sz <- tsne.image.size(facet1=region.count, facet2=gene.count, display.width=display.width)

  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,clusters.selected()$cluster,user.genes(),input$top.N))

  renderCacheImage(tsne.plot, glue("tsne_global_cluster_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)


output$tsne.global.cluster.label.dl <- downloadHandler(filename="tsne.zip", 
                                                       content= function(file) {
                                                         tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=FALSE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.xy(), return.closure = TRUE)()
                                                         send.zip(tsne.plot, 'tsne', file, c('tSNE-plot.R','dv_label.R'))
                                                       })

#########################################################
# show global tSNE with just the filtered subclusters 
output$tsne.global.subcluster.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())

  if (is.filtered()) {
    tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=TRUE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.subcluster.xy())
    
    region.count <- nrow(regions.selected())
    #    gene.or.ic.count <- length(na.omit(c(user.genes(), selected.components()$ic.number)))
    gene.count <- facet2.count()
    #    if (tx.facet2() && gene.or.ic.count==0) { gene.or.ic.count <- length(user.genes()) }
    
    display.width <- img.size.round(session$clientData[[glue("output_tsne.global.subcluster.label_width")]])
    img.sz <- tsne.image.size(facet1=region.count, facet2=gene.count, display.width=display.width)
    
  } else {
    tsne.plot <- function(progress) plot.text("Highlight one or more regions, classes or clusters to begin subcluster analysis")
    h <- img.size.round(session$clientData[[glue("output_tsne.global.subcluster.label_height")]])
    img.sz <- list(height=h,width=h)
  } 
  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,subclusters.selected()$subcluster,user.genes(),input$top.N))
  
  renderCacheImage(tsne.plot, glue("tsne_global_subcluster_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)

output$tsne.global.subcluster.label.dl <- downloadHandler(filename="tsne.zip", 
                                                          content= function(file) {
                                                            tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=TRUE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.subcluster.xy(), return.closure = TRUE)()
                                                            send.zip(tsne.plot, 'tsne', file, c('tSNE-plot.R','dv_label.R'))
                                                          })


#########################################################
output$tsne.local.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())
  
  if (is.filtered()) {
    tsne.plot <- tsne.label(is.global=FALSE, show.subclusters = TRUE, show.cells=TRUE, show.bags=input$opt.show.bags, diff.genes = expr.subcluster.local.xy())
  
    gene.count <- facet2.count()
    #    gene.or.ic.count <- length(na.omit(c(user.genes(), selected.components()$ic.number)))
    cluster.count <- nrow(clusters.selected())
    #    if (tx.facet2() && gene.or.ic.count==0) { gene.or.ic.count <- length(user.genes()) }
    
    display.width <- img.size.round(session$clientData[[glue("output_tsne.local.label_width")]])
    img.sz <- tsne.image.size(facet1=cluster.count, facet2=gene.count, display.width=display.width)

  } else {
    tsne.plot <- function(progress) plot.text("Highlight one or more regions, classes or clusters to begin subcluster analysis")
    h <- img.size.round(session$clientData[[glue("output_tsne.local.label_height")]])
    img.sz <- list(height=h,width=h)
  }
  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,subclusters.selected()$subcluster, user.genes(),input$top.N))
  
  renderCacheImage(tsne.plot, glue("tsne_local_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)

output$tsne.local.label.dl <- downloadHandler(filename="tsne.zip", 
                                              content= function(file) {
                                                tsne.plot <- tsne.label(is.global=FALSE, show.subclusters = TRUE, show.cells=TRUE, show.bags=input$opt.show.bags, diff.genes = expr.subcluster.local.xy(), return.closure = TRUE)()
                                                send.zip(tsne.plot, 'tsne', file, c('tSNE-plot.R', 'dv_label.R'))
                                              })
