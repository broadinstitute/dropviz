source("dv_label.R")

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
tsne.cluster.bag.data <- 
  readRDS(glue("{prep.dir}/tsne/global.clusters.bags.Rdata"))

tsne.subcluster.bag.data <- 
  readRDS(glue("{prep.dir}/tsne/global.subclusters.bags.Rdata"))

tsne.local.bag.data <- 
  readRDS(glue("{prep.dir}/tsne/local.subclusters.bags.Rdata")) 

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
      mutate(cx.gg=factor(ifelse(is.na(cx.selected),NA,as.character(subcluster)), levels=levels(polygon$subcluster))) %>%
      inner_join(cluster.nms, by=c('exp.label','subcluster')) %>%
      dplyr::rename(cx=subcluster, cx.disp=subcluster.disp)
  } else {
    stopifnot(kind=='local')
    major.clusters <- select(clusters, exp.label, cluster) %>% unique
    x <- filter(polygon, exp.label %in% regions$exp.label) %>% inner_join(region.nms, by='exp.label') %>%
      right_join(select(major.clusters, exp.label, cluster), by=c('exp.label','cluster')) %>%
      left_join(select(clusters, exp.label, subcluster) %>% mutate(cx.selected=TRUE), by=c('exp.label','subcluster')) %>%
      mutate(cx.gg=factor(ifelse(is.na(cx.selected),NA,as.character(subcluster)), levels=levels(polygon$subcluster))) %>%
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

global.xy.cluster.selected <- reactive({
log.reactive("fn: global.xy.cluster.selected")
  left_join(select(clusters.selected(),exp.label, cluster), global.xy(), by=c('exp.label','cluster')) %>%
    inner_join(region.names(), by='exp.label') %>%
    dplyr::rename(cx=cluster)
})
global.xy.subcluster.selected <- reactive({
log.reactive("fn: global.xy.subcluster.selected")
  left_join(select(subclusters.selected(),exp.label, subcluster), global.xy(), by=c('exp.label','subcluster')) %>%
    inner_join(region.names(), by='exp.label') %>% 
    select(-cluster) %>% dplyr::rename(cx='subcluster')
})

# returns the xys for the region corresponding with the currently selected cluster
global.xy.current.cluster <- reactive({
log.reactive("fn: global.xy.current.cluster")
  filter(global.xy.cluster.selected(), exp.label==current.cluster()$exp.label & cx==clusters.selected()$cluster[current.cluster()])
}) 
global.xy.current.subcluster <- reactive({
log.reactive("fn: global.xy.current.subcluster")
  filter(global.xy.subcluster.selected(), 
         exp.label==current.subcluster()$exp.label & 
           cx==current.subcluster()$subcluster)
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

# returns the sum of all of the log normal transcript counts for all user.genes
psum.amounts <- function(cx, amounts) {
  amounts$total <- log(rowSums(as_tibble(lapply(as.list(amounts), function(a) exp(a)))))
  cbind(cx, alpha=amounts$total)
}
cluster.transcript.amounts <- reactive({
  cx <- select(clusters.selected(), exp.label, cx=cluster)
  amounts <- select(clusters.selected(), ends_with('.log.target.u'))
  psum.amounts(cx, amounts)
})
subcluster.transcript.amounts <- reactive({
  cx <- select(subclusters.selected(), exp.label, cx=subcluster)
  amounts <- select(subclusters.selected(), ends_with('.log.target.u'))
  psum.amounts(cx, amounts)
})


################################################################
# Draw tSNE

alpha.na2zero <- function(df) mutate(df, alpha=ifelse(is.na(alpha),0,alpha))

# just a short hand to join with alpha only if lhs has data and then replace all NA alphas with zero
left_join_alpha <- function(lhs, rhs, by) {
  if (nrow(lhs) > 0) {
    left_join(lhs, rhs, by=by) %>% alpha.na2zero()
  } else {
    mutate(lhs, alpha=double())
  }
}

# When comps are selected in plot, then data is automatically limited to corresponding cluster
limit.cluster <- function(df, comps) {
  return(df)
  # if (nrow(comps)>0) {
  #   filter(df, cluster==first(comps$cluster) & exp.label==first(comps$exp.label))
  # } else {
  #   df
  # }
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
  c(input$opt.cluster.disp, input$use.common.name, downsample(), input$opt.downsampling.method, input$opt.region.disp, input$opt.plot.label, input$use.bag.plot, input$opt.expr.size)
})

## returns a function that returns either (1) a ggplot object to draw a tsne plot or (2) a closure of the ggplot object and all dependencies
## 
tsne.label <- function(is.global=TRUE, show.subclusters=FALSE, show.cells=TRUE, show.bags=FALSE, diff.genes=tibble(), comps=tibble(), return.closure=FALSE) {
  function(progress=NULL) {
    stopifnot(is.global || show.subclusters) # can't show clusters on local tsne
    
    # for both global.xy (the positions of each cell) and global.[sub]cluster.avg.xy (the center position of each [sub]cluster),
    # limit to selected [sub]clusters and add pretty region name
    write.log(glue("Building tsne.label - is.global={is.global} show.subclusters={show.subclusters} show.bags={show.bags} show.cells={show.cells} diff.genes={nrow(diff.genes)} comps={nrow(comps)}"))
    xy.data <-
      (if (show.cells) {
        if (!is.null(progress)) progress$inc(0.2, detail="Reading XY data")
        if (show.subclusters) {
          if (is.global) {
            global.xy.subcluster.selected()
          } else {
            local.xy.selected()
          }
        } else {
          global.xy.cluster.selected()
        }
      } else {
        tibble()
      }) %>% limit.cluster(comps)

    if (!is.null(progress) && nrow(xy.data)>0) progress$inc(0.2, detail=glue("Read {nrow(xy.data)} cells"))
    
    # labels
    label.data <-
      (
        if (show.subclusters) {
          if (is.global) {
            filter(global.selected.sub.center(), !is.na(cx.gg))
          } else {
            # sometimes there is no local cluster map. If so, then the centers will be (0,0) and there's no points and no bag/loop.
            mutate(filter(local.selected.sub.center(), !is.na(cx.gg)),
                   cx.disp=ifelse(x==0 & y==0, 'No Data', as.character(cx.disp)))
            
          }
        } else {
          na.omit(global.selected.center())
        }
      ) %>% limit.cluster(comps)

    if (show.subclusters) {
      if (is.global) {
        bag.data <- global.selected.sub.bag() %>% limit.cluster(comps)
        loop.data <- global.selected.sub.loop() %>% limit.cluster(comps)
        center.data <- global.selected.sub.center() %>% limit.cluster(comps)
      } else {
        bag.data <- local.selected.sub.bag() %>% limit.cluster(comps)
        loop.data <- local.selected.sub.loop() %>% limit.cluster(comps)
        center.data <- local.selected.sub.center() %>% limit.cluster(comps)
      }
    } else {
      bag.data <- global.selected.bag() %>% limit.cluster(comps)
      loop.data <- global.selected.loop() %>% limit.cluster(comps)
      center.data <- global.selected.center() %>% limit.cluster(comps)
    }
    if (!is.null(progress)) progress$inc(0.2, detail=glue("{nrow(center.data)} cx bag data"))
    
    # alpha is either fixed or set by transcript amounts.
    # if fixed, then alpha range is set by scale further below.
    # if user.genes are specified, then sum all of the amounts per cx and use that as the alpha for colors
    opt.tx.alpha <- !is.null(input$user.genes) && !is.null(input$opt.tx.alpha) && input$opt.tx.alpha
    if (opt.tx.alpha) {
      tx.alpha <- (
        if (show.subclusters) {
          subcluster.transcript.amounts() 
        } else {
          cluster.transcript.amounts()
        }
      ) %>% group_by(exp.label) %>% top_n(as.integer(input$top.N), alpha) 
      
      label.data <- left_join_alpha(label.data, tx.alpha, by=c('exp.label','cx')) 
      center.data <- left_join_alpha(center.data, tx.alpha, by=c('exp.label','cx')) 
      bag.data <- left_join_alpha(bag.data, tx.alpha, by=c('exp.label','cx'))  %>% mutate(alpha=pmax(0,alpha-0.5))
      loop.data <- left_join_alpha(loop.data, tx.alpha, by=c('exp.label','cx'))  %>% mutate(alpha=pmax(0,alpha-1))
      xy.data <- left_join_alpha(xy.data, tx.alpha, by=c('exp.label','cx'))   %>% mutate(alpha=pmax(0,alpha-1))
      
      if (!is.null(progress)) progress$inc(0.2, detail=glue("Computing alpha for {input$user.genes}"))
    }
    
    
    # diff exp & components are row facets
    facet2.vals <- c()
    if (nrow(diff.genes)>0) {
      diff.genes$facet2.gg <- diff.genes$gene
      facet2.vals <- c(facet2.vals, unique(diff.genes$facet2.gg))
    }
    if (nrow(comps)>0) {
      comps$facet2.gg <- comps$IC
      facet2.vals <- c(facet2.vals, unique(comps$facet2.gg))
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
    diff.data <- diff.genes
    comp.data <- comps
    opt.expr.size <- input$opt.expr.size
    opt.show.cells <- show.cells
    opt.show.bags <- show.bags
    opt.global <- is.global
    opt.plot.label <- input$opt.plot.label
    opt.horiz.facet <- nrow(diff.data)>0 || nrow(comp.data)>0

    p.func <- function() {
      require(ggplot2);
      geom_blank_tsne <- geom_blank(data=data.frame(region.disp=character(),facet.gg=character(),facet2.gg=character()))
      tsne.gg <- ggplot() + scale_fill_discrete(guide="none", na.value='lightgray') + scale_size(guide="none", range=c(0.1,opt.expr.size))
      
      tsne.color.scale <- (
        if (nrow(comp.data)>0) {
          scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint=0, limits=c(-max(comp.data$weight),max(comp.data$weight)))
        } else {
          scale_color_discrete(guide="none", na.value='lightgray') 
        }
      )

      xy.gg <- (
        if (opt.show.cells) {
          if (nrow(comp.data)>0) {
            geom_point(data=xy.data, aes(x=V1,y=V2), color='grey', alpha=0.25, size=xy.cell.size()) 
          } else {
            if (opt.tx.alpha) {
              geom_point(data=xy.data, aes(x=V1,y=V2,color=cx, alpha=alpha), size=xy.cell.size()) 
            } else {
              geom_point(data=xy.data, aes(x=V1,y=V2,color=cx), alpha=0.25, size=xy.cell.size()) 
            }
          }
        } else {
          geom_blank_tsne
        }
      )

      label.gg <- (
        if (nrow(comp.data)>0) {
          geom_label(data=label.data, aes(x=x,y=y,label=as.character(cx.disp)), color='grey')
        } else {
          if (opt.tx.alpha) {
            dv_label(data=label.data, aes(x=x,y=y,color=cx.gg,label=as.character(cx.disp),alpha=alpha))
          } else {
            dv_label(data=label.data, aes(x=x,y=y,color=cx.gg,label=as.character(cx.disp)))
          }
        }
      )
      
      diff.gg <- (
        if (nrow(diff.data)>0) {
          geom_point(data=diff.data, aes(x=V1, y=V2, size=transcripts), color='black', alpha=0.2)
        } else {
          geom_blank_tsne
        }
      )

      comp.gg <- (
        if (nrow(comp.data)>0) {
          geom_point(data=comp.data, aes(x=V1, y=V2, color=weight), size=0.75, alpha=1)
        } else {
          geom_blank_tsne
        }
      )      
      
      if (opt.show.bags & nrow(comp.data)==0) {
        if (opt.tx.alpha) {
          alpha.range <- scale_alpha_continuous(guide="none",range=c(0,1), trans=scales::trans_new("sqr", function(x) x^2, function(x) sqrt(x)))
          bag.gg <- geom_polygon(data=filter(bag.data, !is.na(cx.gg)), aes(x=x,y=y,fill=cx.gg,group=cx, alpha=alpha))
          loop.gg <- geom_polygon(data=filter(loop.data, !is.na(cx.gg)), aes(x=x,y=y,fill=cx.gg,group=cx, alpha=alpha))
          center.gg <- geom_point(data=filter(center.data, !is.na(cx.gg)), aes(x=x,y=y,color=cx.gg, fill=cx.gg,alpha=alpha), size=3)
        } else {
          bag.gg <- geom_polygon(data=bag.data, aes(x=x,y=y,fill=cx.gg,group=cx), alpha=0.4)
          loop.gg <- geom_polygon(data=loop.data, aes(x=x,y=y,fill=cx.gg,group=cx), alpha=0.2)
          center.gg <- geom_point(data=center.data, aes(x=x,y=y,color=cx.gg, fill=cx.gg), size=3)
          alpha.range <- scale_alpha()
        }
      } else {
        bag.gg <- geom_polygon(data=bag.data, aes(x=x,y=y,group=cx), fill='grey', alpha=0.2) 
        loop.gg <- geom_polygon(data=loop.data, aes(x=x,y=y,group=cx), fill='grey', alpha=0.1) 
        center.gg <- geom_blank_tsne
        alpha.range <- scale_alpha()
      }
      
      facet.label.gg <- (
        if (opt.horiz.facet) {
          geom_text(data=facet.label.data, aes(x=x, y=y, label=facet2.gg), hjust="left")
        } else {
          geom_blank_tsne
        }
      )

      p <- tsne.gg + tsne.color.scale + xy.gg + loop.gg + bag.gg + alpha.range + center.gg + diff.gg + comp.gg + label.gg + facet.label.gg 
      
      if (opt.global) {
        if (opt.horiz.facet) {
          plot.gg <- p + facet_grid(facet2.gg~region.disp)
        } else {
          plot.gg <- p + facet_wrap(~region.disp, ncol=4)
        }
      } else {
        if (opt.horiz.facet) {
          plot.gg <- p + facet_grid(facet2.gg~region.disp+facet.gg)
        } else {
          plot.gg <- p + facet_wrap(~region.disp+facet.gg, ncol=ifelse(opt.global,3,4))
        }
      }
      plot.gg
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
  isTruthy(input$tissue) || isTruthy(input$cell.class) || isTruthy(input$cell.cluster) || isTruthy(input$cell.type)
}

# Set the image size to create square facets. The display is divided into at most 4 facets across when wrapping (e.g. 250x250 for a 1000px region).
# If there's a second facet, then a facet_grid display is used only if the 2nd facet is less than MAX_REGIONS (otherwise each facet would be too small)
tsne.image.size <- function(facet1, facet2, display.width) {
  if (facet2 > 0 && facet1 <= MAX_REGIONS) {
    # grid: 
    facet.wide = (if (facet1==1) display.width %/% 2 else display.width)
    facet.high <- if (facet2==1) display.width %/% 2 else facet2 * (display.width %/% 4)
  } else {
    # wrap
    facet.wide.count <- min(4, facet1)
    facet.wide <- ifelse(facet.wide.count==1, display.width %/% 2, display.width)
    facet.high <- (facet.wide %/% facet.wide.count) * (((facet1-1) %/% facet.wide.count)+1)
  }
  
  return(list(width=facet.wide, height=facet.high))
}

#########################################################
# show global tSNE with just the filtered clusters 
output$tsne.global.cluster.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())
  
  tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=FALSE, show.cells=downsample()>0, show.bags = input$opt.show.bags, diff.genes=expr.xy())
  
  region.count <- nrow(regions.selected())
  gene.count <- nrow(cluster.markers.selected())

  img.sz <- tsne.image.size(facet1=region.count, facet2=gene.count, display.width=img.size.round(session$clientData[[glue("output_tsne.global.cluster.label_width")]]))

  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,clusters.selected()$cluster,cluster.markers.selected()$gene,input$user.genes,input$top.N,input$opt.show.bags,input$opt.tx.alpha))

  renderCacheImage(tsne.plot, glue("tsne_global_cluster_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)


output$tsne.global.cluster.label.dl <- downloadHandler(filename="tsne.zip", 
                                                       content= function(file) {
                                                         tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=FALSE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.xy(), return.closure = TRUE)()
                                                         send.zip(tsne.plot, 'tsne', file)
                                                       })

#########################################################
# show global tSNE with just the filtered subclusters 
output$tsne.global.subcluster.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())

  if (is.filtered()) {
    tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=TRUE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.subcluster.xy())
    
    region.count <- nrow(regions.selected())
    gene.or.ic.count <- length(na.omit(c(subcluster.markers.selected()$gene, selected.components()$ic.number)))
    
    img.sz <- tsne.image.size(facet1=region.count, facet2=gene.or.ic.count, display.width=img.size.round(session$clientData[[glue("output_tsne.global.subcluster.label_width")]]))
    
  } else {
    tsne.plot <- function(progress) plot.text("Highlight one or more regions, classes or clusters to begin subcluster analysis")
    h <- img.size.round(session$clientData[[glue("output_tsne.global.subcluster.label_height")]])
    img.sz <- list(height=h,width=h)
  } 
  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,subclusters.selected()$subcluster,subcluster.markers.selected()$gene,input$user.genes,input$top.N,input$opt.show.bags,input$opt.tx.alpha))
  
  renderCacheImage(tsne.plot, glue("tsne_global_subcluster_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)

output$tsne.global.subcluster.label.dl <- downloadHandler(filename="tsne.zip", 
                                                          content= function(file) {
                                                            tsne.plot <- tsne.label(is.global=TRUE, show.subclusters=TRUE, show.cells=(downsample()>0), show.bags = input$opt.show.bags, diff.genes=expr.subcluster.xy(), return.closure = TRUE)()
                                                            send.zip(tsne.plot, 'tsne', file)
                                                          })


#########################################################
output$tsne.local.label <- renderImage({
  progress <- shiny.progress('t-SNE')
  if (!is.null(progress)) on.exit(progress$close())
  
  if (is.filtered()) {
    tsne.plot <- tsne.label(is.global=FALSE, show.subclusters = TRUE, show.cells=TRUE, show.bags=input$opt.show.bags, diff.genes = expr.subcluster.local.xy(), comps=selected.component.cell.weights.xy())
  
    gene.or.ic.count <- length(na.omit(c(subcluster.markers.selected()$gene, selected.components()$ic.number)))
    cluster.count <- nrow(clusters.selected())
    
    img.sz <- tsne.image.size(facet1=cluster.count, facet2=gene.or.ic.count, display.width=img.size.round(session$clientData[[glue("output_tsne.local.label_width")]]))

  } else {
    tsne.plot <- function(progress) plot.text("Highlight one or more regions, classes or clusters to begin subcluster analysis")
    h <- img.size.round(session$clientData[[glue("output_tsne.local.label_height")]])
    img.sz <- list(height=h,width=h)
  }
  key.str <- digest(c(tsne.disp.opts(),regions.selected()$exp.label,subclusters.selected()$subcluster, subcluster.markers.selected()$gene, selected.components()$ic.number,input$user.genes,input$top.N,input$opt.show.bags,input$opt.tx.alpha))
  
  renderCacheImage(tsne.plot, glue("tsne_local_label_{key.str}"), img.sz$width, img.sz$height, progress=progress)
}, deleteFile = FALSE)

output$tsne.local.label.dl <- downloadHandler(filename="tsne.zip", 
                                              content= function(file) {
                                                tsne.plot <- tsne.label(is.global=FALSE, show.subclusters = TRUE, show.cells=TRUE, show.bags=input$opt.show.bags, diff.genes = expr.subcluster.local.xy(), comps=selected.component.cell.weights.xy(), return.closure = TRUE)()
                                                send.zip(tsne.plot, 'tsne', file)
                                              })
