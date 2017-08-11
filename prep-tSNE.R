# pre-generate tSNEs data for each cluster, sub-cluster

source("global.R")
tsne.dir <- glue("{prep.dir}/tsne")

# reads the XY tSNE coordinates for each cell. Returns as a tibble.
read.tSNE.xy <- function(fname) {
  tSNE.df <- as.data.frame(readRDS(fname)) # Must first read as data.frame to get rownames
  tSNE.xy <- as_tibble(tSNE.df)
  tSNE.xy$cell <- rownames(tSNE.df)
  tSNE.xy
}

# print the plot p to the file fn at different sizes and resolutions
#plot.dims <- c(500, 1500)
# plot.multisize <- function(fn, p) {
#   sapply(plot.dims, function(sz) {
#     fn <- paste0(fn,"_",sz,".png")
#     log("Plotting: ", fn)
#     png(fn, width=sz,height=sz,units="px",res=sz*.15)
#     print(p)
#     dev.off()
#   })
#   
# }

dlply(experiments, .(exp.label), function(exp) {
  log("Generating tSNE data for  ",exp$exp.label)
  out.dir <- sprintf("%s/%s",tsne.dir, exp$exp.label)
  dir.create(out.dir, recursive = TRUE, showWarnings=FALSE)
  
  # GLOBAL tSNE
  tSNE.xy <- read.tSNE.xy(sprintf("%s/tSNE/%s_tSNExy.RDS", exp$exp.dir,exp$base))

  # convert factor to table
  cluster.cell.assign <- readRDS(sprintf("%s/assign/%s.cluster.assign.RDS", exp$exp.dir, exp$base))
  cluster.cell.assign.tbl <- tibble(cell=names(cluster.cell.assign), cluster=factor(cluster.cell.assign, levels=levels(cell.types$cluster)))
  
  subcluster.cell.assign <- readRDS(sprintf("%s/assign/%s.subcluster.assign.RDS", exp$exp.dir, exp$base))
  subcluster.cell.assign.tbl <- tibble(cell=names(subcluster.cell.assign), subcluster=factor(subcluster.cell.assign, levels=levels(cell.types$subcluster)))
  
  # merge cluster and subcluster into a single table
  # cell subcluster cluster
  # <chr>     <fctr>  <fctr>
  # 1  P60STRRep1P1_AGCCGCTTAATA        2-3       2
  # 2  P60STRRep1P1_GTGTCGTCCGCT        2-3       2
  # 3  P60STRRep1P1_ACTCTACCAAAT        2-5       2
  # 4  P60STRRep1P1_CGGTGTGACTAC        2-5       2
  cell.assign <- full_join(subcluster.cell.assign.tbl, cluster.cell.assign.tbl, by='cell')
  
  tSNE <- inner_join(tSNE.xy, cell.assign, by='cell')
  
  global.xy.fn <- sprintf("%s/global.xy.RDS", out.dir)
  log("Writing ", global.xy.fn)
  saveRDS(tSNE, global.xy.fn)

  # # generate different size images
  # lapply(c(500, 1500), function(sz) {
  #   log("Creating ", sz, "x", sz, " images")
  #   
  #   # global tSNE of major cluster labels
  #   plot.multisize(sprintf("%s/%s_class", out.dir, exp$exp.label), 
  #                  ggplot(tSNE, aes(x=V1, y=V2, color=cluster)) + geom_point(alpha=0.25, size=0.5) + ggtitle(exp$exp.label))
  # 
  #   # global tSNE of subcluster labels
  #   plot.multisize(sprintf("%s/%s_type", out.dir, exp$exp.label),
  #                  ggplot(filter(tSNE,!is.na(subcluster)), aes(x=V1, y=V2, color=subcluster)) + geom_point(alpha=0.25, size=0.5))

  # for each major cluster, generate tSNE with subcluster labels
  lapply(unique(filter(cell.types, exp.label==exp$exp.label)$cluster), function(cn) {
      log("Processing cluster ",cn)
      fn <- sprintf("%s/cluster%s/%s.cluster%s.auto.tSNExy.RDS", exp$exp.dir,cn,exp$base,cn)
      if (file.exists(fn)) {
        
        # local coords of cell
        tSNE.xy <- read.tSNE.xy(fn)
        
        # subcluster labeling of each cell
        # V1         V2                      cell subcluster cluster
        # <dbl>      <dbl>                     <chr>     <fctr>  <fctr>
        # 1  -24.879830 -6.0923347 P60STRRep1P1_AGTGCAAACTGT         NA       8
        # 2   11.933133 -2.0367023 P60STRRep1P1_AGTCATTTCATA        8-1       8
        # 3    5.898253  0.4181663 P60STRRep1P1_GTCCTTCCTGGG        8-1       8
        tSNE.local <- inner_join(tSNE.xy, filter(cell.assign, cluster==cn & !is.na(subcluster)), by='cell')
        
        local.xy.fn <- sprintf("%s/cluster%s.xy.RDS", out.dir, cn)
        log("Writing ", local.xy.fn)
        saveRDS(tSNE.local, local.xy.fn)

        ## if (nrow(tSNE.local) > 0) {
          
        ##   plot.multisize(sprintf("%s/%s_%s", out.dir, exp$exp.label, cn),
        ##                  ggplot(tSNE.local, aes(x=V1, y=V2, color=subcluster)) + geom_point(size=1) + ggtitle(sprintf("Cluster %s",cn))
        ##   ) # + theme(legend.position = 'bottom') # + facet_wrap(~cluster)

        ##   expr.fn <- sprintf("%s/cluster%s/%s_cluster%s.DGE.RDS", exp$exp.dir,cn,exp$base,cn)
          
        ##   if (file.exists(expr.fn)) {
        ##     expr <- readRDS(expr.fn) # rows genes, cols cells in this major cluster
        ##     lapply(all.genes, function(gene.name) {

        ##       if (gene.name %in% rownames(expr)) {
        ##         gene.expr <- expr[gene.name,]
                
        ##         # #                         cell transcript.count
        ##         # <chr>            <dbl>
        ##         # 1  P60STRRep1P1_AGTGCAAACTGT                0
        ##         # 2  P60STRRep1P1_AGTCATTTCATA                2
        ##         # 3  P60STRRep1P1_GTCCTTCCTGGG                0
        ##         # 4  P60STRRep1P1_GTGGATTTTCCT                0
        ##         gene.expr.tbl <- tibble(cell=names(gene.expr), transcript.count=gene.expr)
                
        ##         # V1         V2                      cell subcluster cluster transcript.count
        ##         # <dbl>      <dbl>                     <chr>     <fctr>  <fctr>            <dbl>
        ##         # 1  -24.879830 -6.0923347 P60STRRep1P1_AGTGCAAACTGT         NA       8                0
        ##         # 2   11.933133 -2.0367023 P60STRRep1P1_AGTCATTTCATA        8-1       8                2
        ##         # 3    5.898253  0.4181663 P60STRRep1P1_GTCCTTCCTGGG        8-1       8                0
        ##         # 4  -24.639872 -5.4817099 P60STRRep1P1_GTGGATTTTCCT         NA       8                0
        ##         tSNE.local.expr <- left_join(filter(tSNE.local, !is.na(subcluster)), gene.expr.tbl, by='cell') %>%
        ##           mutate(is_expressed=transcript.count>0)
                
        ##         plot.multisize(sprintf("%s/%s_%s_%s", out.dir, exp$exp.label, cn, gene.name),
        ##                        ggplot(tSNE.local.expr, aes(x=V1, y=V2, size=transcript.count, color=subcluster, shape=is_expressed)) + geom_point() + ggtitle(sprintf("Cluster %s / Gene %s", cn, gene.name))
        ##         )
        ##       } else {
        ##         log(gene.name," is missing from ", expr.fn, warn = TRUE)
        ##       }
        ##       gene.name
        ##     })
        ##   } else {
        ##     log(expr.fn," not found", warn=TRUE)
        ##   }
        ## } else {
        ##   log("No subcluster cell assignments for cluster ", cn, warn=TRUE)
        ## }
      
      } else {
        log(fn," not found", warn=TRUE)
      }
      
    })
  # })
  
})


## generate bag plot data
## Custom bagplot hack
##
## Using code from aplpack https://cran.r-project.org/web/packages/aplpack/index.html
## that was ripped out and made into a geom_bag for ggplot https://gist.github.com/benmarwick/00772ccea2dd0b0f1745 (which I don't use)
## this code stores the bag plotting information (polygons) for future display

library(tibble)
library(glue)
library(plyr)
library(dplyr)
library(ggplot2)
source("utils/bag_functions.r")

source("global.R")

mk.bag.xy <- function(xy, grouping, this.exp.label) {
  dlply(xy, grouping, function(df) {
    this.group <- first(df[[grouping]])
    title <- paste(this.exp.label,this.group)
    cat(title,"\n")
    
    ## for some reason, there are clusters referenced in the global coordinates that do not have
    ## cluster subdirs - probably because most are too small.
    if ((grouping == 'cluster' && nrow(filter(cluster.names_, exp.label==this.exp.label & cluster==this.group))>0) ||
        (grouping == 'subcluster' && nrow(filter(subcluster.names_, exp.label==this.exp.label & subcluster==this.group))>0)) {
      m <- data.matrix(df[,c('V1','V2')])
      
      hulls <- hulls_for_bag_and_loop(m)
      
      # close the loop by repeating the first coord at the end
      hulls_closed <- lapply(hulls[1:2], function(i) data.frame(rbind(i, i[1, ] ), row.names = NULL ))
      
      the_loop <-  setNames(data.frame(hulls_closed$hull.loop), nm = c("x", "y"))
      the_bag <-   setNames(data.frame(plothulls_(m, fraction = 0.5)), nm = c("x", "y"))
      
      # get the center, which are new coords not in the original dataset
      center <-  setNames(data.frame(matrix(hulls$center, nrow = 1)), nm = c("x", "y"))
      
      list(loop=the_loop, bag=the_bag, center=center)
    } else {
      warning("Skipping ",title,", points=",nrow(df))
      list()
    }
    
  })
  
}

mk.subc.bag.data <- function() {
  bag.data <-
    dlply(cluster.names_, .(exp.label, cluster), function(df) {
      fn <- glue("{tsne.dir}/{df$exp.label}/cluster{df$cluster}.xy.RDS")
      if (file.exists(fn)) {
        xy <- readRDS(fn)
        mk.bag.xy(xy, 'subcluster', df$exp.label)
      } else {
        warning("Skipping missing file ", fn)
        list()
      }
    })

  list(
    loops=ldply(bag.data, function(exp.cl) ldply(exp.cl, function(x) x$loop)) %>% as_tibble,
    bags=ldply(bag.data, function(exp.cl) ldply(exp.cl, function(x) x$bag)) %>% as_tibble,
    centers=ldply(bag.data, function(exp.cl) ldply(exp.cl, function(x) x$center)) %>% as_tibble)
}

mk.bag.data <- function(kind) {
  bag.data <- 
    dlply(experiments, .(exp.label), function(exp) {
      xy <- readRDS(glue("{tsne.dir}/{exp$exp.label}/global.xy.RDS"))
      mk.bag.xy(xy, kind, exp$exp.label)
    })
  
  list(
    loops=ldply(bag.data, function(exp) ldply(exp, function(x) x$loop)) %>% as_tibble,
    bags=ldply(bag.data, function(exp) ldply(exp, function(x) x$bag)) %>% as_tibble,
    centers=ldply(bag.data, function(exp) ldply(exp, function(x) x$center)) %>% as_tibble)
}  
  
bagdata <- mk.bag.data('cluster')
saveRDS(bagdata, file=glue("{tsne.dir}/global.clusters.bags.Rdata"))

if (FALSE) {
  bagdata <- readRDS(glue("{tsne.dir}/global.clusters.bags.Rdata"))
  ggplot() +
    geom_polygon(data=bagdata$loops, mapping=aes(x=x,y=y, fill=cluster), alpha=0.2) +
    geom_polygon(data=bagdata$bags, mapping=aes(x=x,y=y, fill=cluster), alpha=0.3) +
    geom_point(data=bagdata$centers, mapping=aes(x=x,y=y, color=cluster), size=3) +
    facet_wrap(~exp.label) + scale_fill_discrete(guide="none") + scale_color_discrete(guide="none")
}

bagdata <- mk.bag.data('subcluster')
saveRDS(bagdata, file=glue("{tsne.dir}/global.subclusters.bags.Rdata"))

if (FALSE) {
  bagdata <- readRDS(glue("{tsne.dir}/global.subclusters.bags.Rdata"))
  ggplot() +
    geom_polygon(data=bagdata$loops, mapping=aes(x=x,y=y, fill=subcluster), alpha=0.2) +
    geom_polygon(data=bagdata$bags, mapping=aes(x=x,y=y, fill=subcluster), alpha=0.3) +
    geom_point(data=bagdata$centers, mapping=aes(x=x,y=y, color=subcluster)) +
    facet_wrap(~exp.label) + scale_fill_discrete(guide="none") + scale_color_discrete(guide="none")
}


bagdata <- mk.subc.bag.data()


# some sub-classes have no data. Create an indicator of a center of (0,0)
# the chance that a true center is exactly (0,0) is probably infinitesimally small
add.zero <- function(df) {
  right_join(df, select(cell.types, -class), by=c('exp.label','cluster','subcluster')) %>% 
    mutate(x=ifelse(is.na(x),0,x),y=ifelse(is.na(y),0,y))
}

bagdata$centers <- add.zero(bagdata$centers)
bagdata$loops <- add.zero(bagdata$loops)
bagdata$bags <- add.zero(bagdata$bags)

saveRDS(bagdata, file=glue("{tsne.dir}/local.subclusters.bags.Rdata"))

if (FALSE) {
  bagdata <- readRDS(glue("{tsne.dir}/local.subclusters.bags.Rdata"))
  ggplot() +
    geom_polygon(data=bagdata$loops, mapping=aes(x=x,y=y, fill=subcluster), alpha=0.2) +
    geom_polygon(data=bagdata$bags, mapping=aes(x=x,y=y, fill=subcluster), alpha=0.3) +
    geom_point(data=bagdata$centers, mapping=aes(x=x,y=y, color=subcluster)) +
    facet_wrap(~exp.label+cluster) + scale_fill_discrete(guide="none") + scale_color_discrete(guide="none")
}

