    # DropViz - rank plot
    #
    #   This R code along with the .Rdata bundled in the zip file is
    #   used to generate a rank plot of the N (sub)clusters that have
    #   the highest expression for the genes of interest.
    #
    #   To execute the code, (1) set your working directory to the
    #   directory containing this file, (2) load the data with
    #   load("rank.Rdata") or you may be able to double-click the file
    #   to load the data into your working environment, (3) source
    #   this file, i.e. source('rank-plot.R') (4) print, display or save
    #   the ggplot object called plot.gg, e.g. print(plot.gg) should
    #   display the scatter plot on your current graphics device.

    # Install from CRAN
    require(ggplot2)
    require(glue)
    require(ggthemes)


    plot.gg <- ggplot(clusters.top, 
                      aes(x=target.sum.per.100k, xmin=target.sum.L.per.100k, 
                          xmax=target.sum.R.per.100k, y=cx.disp, yend=cx.disp)) + 
      geom_point(size=3) + 
      geom_segment(aes(x=target.sum.L.per.100k,xend=target.sum.R.per.100k)) +
      ggtitle(glue("{Kind}s With Highest Expression of {paste0(genes,collapse=' & ')} (top {length(unique(clusters.top$cx.disp))} results)")) + 
      xlab(glue("Transcripts Per 100,000 in {Kind}\n\nThe reported confidence intervals reflect statistical sampling noise (calculated from the binomial distribution,\nand reflecting total number of UMIs ascertained by cluster) rather than cell-to-cell heterogeneity within a cluster")) + 
      ylab("") + 
      rank.facet_grid + 
      xlim(0, max(clusters.top$target.sum.R.per.100k)) +
      theme_few() + 
      theme(plot.title=element_text(size=20,face="bold",hjust=0.5), strip.text=element_text(size=16), axis.text.y=element_text(size=16))
