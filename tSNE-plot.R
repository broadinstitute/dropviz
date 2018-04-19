
      # DropViz - tSNE plot
      #
      #   This R code along with the .Rdata bundled in the zip file is
      #   used to generate the tSNE plots. There are many options, but
      #   the parameters in the .Rdata file correspond to those set a
      #   the time of download. 
      #
      #   To execute the code, (1) set your working directory to the
      #   directory containing this file, (2) load the data with
      #   load("tsne.Rdata") or you may be able to double-click the
      #   file to load the data into your working environment, (3)
      #   source this file, i.e. source('tsne-plot.R') (4) print, display
      #   or save the ggplot object called plot.gg,
      #   e.g. print(plot.gg) should display the tSNE plot on your
      #   current graphics device.
      #
      #   If you are familiar with ggplot then it should be fairly
      #   straightforward to tweak the plot to your preferences. If
      #   you have questions, email dkulp@broadinstitute.org.
      require(ggplot2)                             # install from CRAN
      require(ggthemes)                            # install from CRAN
      if (!exists('dv_label')) source("dv_label.R") # custom ggplot include in zip

      # randomly, but reproducibly, permute the levels of fctr.
      # fctr is either a cluster or subcluster. Without permutation
      # then ggplot uses similar colors for subclusters that are
      # adjacent in order, e.g. 4-1, 4-2, etc.
      cx.permute <- function(fctr) {
        seed.save <- if (exists('.Random.seed')) .Random.seed else NULL
        set.seed(10)
        fctr <- factor(fctr, levels=sample(levels(fctr)))
        .Random.seed <- seed.save
        fctr
      }

      if (opt.tx.heat) {
        scale_color <- function(...) scale_color_brewer(..., palette='YlOrRd')
        scale_fill <- function(...) scale_fill_brewer(..., palette='YlOrRd')
      } else {
        if (length(unique(center.data$cx.gg)) <= 20) {
          # gdocs has a max of 20 discrete colors
          scale_fill <- scale_fill_gdocs
          scale_color <- scale_color_gdocs
        } else {
          scale_fill <- scale_fill_discrete
          scale_color <- scale_color_discrete
        }
      }

      geom_blank_tsne <- geom_blank(data=data.frame(region.disp=character(),facet.gg=character(),facet2.gg=character()))

      tsne.fill.scale <- (
        if (opt.tx.heat) {
          scale_fill(guide="none", na.value='lightgray', limits=levels(center.data$heat))
        } else {
          scale_fill(guide="none", na.value='lightgray')
        }
      )
      
      tsne.color.scale <- (
        if (opt.tx.heat) {
          scale_color(guide=opt.tx.legend, name="Log Expression", na.value='lightgray', limits=levels(center.data$heat))
        } else {
          scale_color(guide="none", na.value='lightgray') 
        }
      )

      # xy.data is (possibly sub-sampled) points of all cells
      xy.gg <- (
        if (opt.show.cells) {
          if (opt.tx.alpha) {
            geom_point(data=xy.data, aes(x=V1,y=V2,color=cx.permute(cx), alpha=alpha), size=opt.xy.cell.size) 
          } else if (opt.tx.heat) {
            geom_point(data=xy.data, aes(x=V1,y=V2,color=heat), alpha=0.25, size=CELL.MIN.SIZE)
          } else {
            geom_point(data=xy.data, aes(x=V1,y=V2,color=cx.permute(cx)), alpha=0.25, size=opt.xy.cell.size) 
          }
        } else {
          geom_blank_tsne
        }
      )

      # label.data are the names of (sub)clusters, positioned at the centroid of each (sub)clusters
      label.gg <- (
        if (opt.tx.heat) {
          dv_label(data=label.data, aes(x=x,y=y,color=heat,label=as.character(cx.disp)), show.legend=FALSE)
        } else {
          dv_label(data=label.data, aes(x=x,y=y,color=cx.permute(cx.gg),label=as.character(cx.disp)), show.legend=FALSE)
        }
      )
      
      # diff.data is expression data on specific genes, either entered by user or selected rows from diff exp table
      diff.gg <- (
        if (nrow(diff.data)>0 && opt.tx.cells) {
          if (opt.cell.display.type=='size') {
            geom_point(data=diff.data, aes(x=V1, y=V2, size=transcripts), color='black', alpha=0.2)
          } else {
            # absent/present
            geom_point(data=diff.data, aes(x=V1, y=V2), size=CELL.MIN.SIZE, color='black', alpha=0.2)
          }
        } else {
          geom_blank_tsne
        }
      )

      # bag.data, loop.data and center.data are the polygon and point data associated with each (sub)cluster
      if (opt.show.bags) {
        if (opt.tx.alpha) {
          alpha.limits <- if (opt.tx.scale=='fixed') c(0,7) else c(0,max(center.data$alpha))
          alpha.range <- scale_alpha_continuous(guide="none", range=c(0,1), limit=alpha.limits)
          bag.gg <- geom_polygon(data=filter(bag.data, !is.na(cx.gg)), aes(x=x,y=y,fill=cx.permute(cx.gg), group=cx, alpha=alpha))
          loop.gg <- geom_polygon(data=filter(loop.data, !is.na(cx.gg)), aes(x=x,y=y,fill=cx.permute(cx.gg), group=cx, alpha=alpha))
          center.gg <- geom_point(data=filter(center.data, !is.na(cx.gg)), aes(x=x,y=y, color=cx.permute(cx.gg), alpha=alpha), size=3)
        } else if (opt.tx.heat) {
          bag.gg <- geom_polygon(data=bag.data, aes(x=x,y=y,fill=heat,group=cx), alpha=0.6)
          loop.gg <- geom_polygon(data=loop.data, aes(x=x,y=y,fill=heat,group=cx), alpha=0.2)
          center.gg <- geom_point(data=center.data, aes(x=x,y=y,color=heat), size=3)
          alpha.range <- scale_alpha()
        } else {
          bag.gg <- geom_polygon(data=bag.data, aes(x=x,y=y,fill=cx.permute(cx.gg),group=cx), alpha=0.4)
          loop.gg <- geom_polygon(data=loop.data, aes(x=x,y=y,fill=cx.permute(cx.gg),group=cx), alpha=0.2)
          center.gg <- geom_point(data=center.data, aes(x=x,y=y,color=cx.permute(cx.gg)), size=3)
          alpha.range <- scale_alpha()
        }
      } else {
        bag.gg <- geom_polygon(data=bag.data, aes(x=x,y=y,group=cx), fill='grey', alpha=0.2) 
        loop.gg <- geom_polygon(data=loop.data, aes(x=x,y=y,group=cx), fill='grey', alpha=0.1) 
        center.gg <- geom_blank_tsne
        if (opt.tx.alpha) {
          alpha.limits <- if (opt.tx.scale=='fixed') c(0,7) else c(0,max(center.data$alpha))
          alpha.range <- scale_alpha_continuous(guide="none", range=c(0,1), limit=alpha.limits)
        } else {
          alpha.range <- scale_alpha()
        }
      }

      # the facet.label.gg is a text label in the top left of each faceted plot.
      # Usually this is the name of the selected gene, a repeat of the facet label on the far right.
      # If there is no horizontal faceting because expression levels of selected genes are summed, then
      # the label is the name of all the genes.
      facet.label.gg <- (
        if (opt.horiz.facet) {
          geom_text(data=facet.label.data, aes(x=x, y=y, label=facet2.gg), hjust="left", show.legend = FALSE)
        } else if (opt.tx.alpha || opt.tx.heat) {
          geom_text(data=tibble(x=min(loop.data$x),
                                y=max(loop.data$y),
                                gene=paste(user.genes(), collapse='+')), aes(x=x,y=y,label=gene), hjust="left", show.legend = FALSE)
        } else {
          geom_blank_tsne
        }
      )
      
      tsne.lim <- 20  # at minimum show box [-20..20]
      xy.limits.gg <- (
        if (opt.show.bags) {
          coord_cartesian(xlim=c(min(c(-tsne.lim, loop.data$x)), max(c(tsne.lim, loop.data$x))), 
                          ylim=c(min(c(-tsne.lim, loop.data$y)), max(c(tsne.lim, loop.data$y))))
        } else { 
          # opt.show.cells
          coord_cartesian(xlim=c(min(c(-tsne.lim, xy.data$V1), na.rm=TRUE), max(c(tsne.lim, xy.data$V1), na.rm=TRUE)), 
                          ylim=c(min(c(-tsne.lim, xy.data$V2), na.rm=TRUE), max(c(tsne.lim, xy.data$V2), na.rm=TRUE)))
        }
      )

      tsne.gg <- ggplot() + tsne.color.scale + tsne.fill.scale + scale_size(guide="none", range=c(0.1,opt.expr.size)) + theme_few() + theme(strip.text.x=element_text(size=20), strip.text.y=element_text(size=14))
      p <- tsne.gg + center.gg + xy.gg + loop.gg + bag.gg + alpha.range + diff.gg + label.gg + facet.label.gg + xy.limits.gg + xlab("V1") + ylab("V2")

      if (opt.global) {
        if (opt.horiz.facet && length(facet2.vals)>1) {
          plot.gg <- p + facet_grid(facet2.gg~region.disp)
        } else {
          plot.gg <- p + facet_wrap(~region.disp, ncol=3)
        }
      } else {
        if (opt.horiz.facet && length(facet2.vals)>1) {
          plot.gg <- p + facet_grid(facet2.gg~region.disp+facet.gg)
        } else {
          plot.gg <- p + facet_wrap(~region.disp+facet.gg, ncol=3)
        }
      }
      plot.gg

