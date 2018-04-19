      # DropViz - scatter plot
      #
      #   This R code along with the .Rdata bundled in the zip file is
      #   used to generate the scatter plots.
      #
      #   To execute the code, (1) set your working directory to the
      #   directory containing this file, (2) load the data with
      #   load("scatter.Rdata") or you may be able to double-click the file
      #   to load the data into your working environment, (3) source
      #   this file, i.e. source('scatter-plot.R') (4) print, display or save
      #   the ggplot object called plot.gg, e.g. print(plot.gg) should
      #   display the scatter plot on your current graphics device.
      #
      #   The primary object used for the plot is called cell.pairs.
      #   It is a tibble containing metacell gene expression in the
      #   target and comparison (sub)cluster(s) for each gene. There
      #   are also columns indicating whether a gene passed criteria,
      #   was selected by the user, etc.  Other parameters, prefixed
      #   with "opt" represent the display settings at the time of
      #   download.
      #
      #   If you are familiar with ggplot then it should be fairly
      #   straightforward to tweak the plot to your preferences. If
      #   you have questions, email dkulp@broadinstitute.org.

      require(ggplot2); require(dplyr)
      
      gene.labels <- (
        if (opt.scatter.gene.labels && any(!is.na(cell.pairs$gene.label))) {
          geom_label(aes(x=log.target.u+.3, label=gene.label), size=3, show.legend = FALSE)
        } else {
          geom_blank()
        }
      )
      
      if (opt.expr.filter %in% c('both','fc','either')) {
        if (opt.expr.filter == 'fc') {
          fc.line1 <- geom_abline(intercept=log(opt.fold.change), slope=1, alpha=0.5, color='grey') 
        } else {
          fc.line1 <- geom_blank()
        }
        fc.line2 <- geom_abline(intercept=-log(opt.fold.change), slope=1, alpha=0.5, color='grey')
      } else {
        fc.line1 <- fc.line2 <- geom_blank()
      }
      
      if (opt.expr.filter %in% c('both','amt','either')) {
        amt.line1 <- geom_hline(yintercept = opt.max.amt.without, color='grey') 
        amt.line2 <- geom_vline(xintercept = opt.min.amt.within, color='grey')
      } else {
        amt.line1 <- amt.line2 <- geom_blank()
      }

      plot.gg <- ggplot(cell.pairs, aes(x=log.target.u, y=log.comparison.u, color=pval.disp, size=size)) + geom_point() + 
        fc.line1 + fc.line2 +
        amt.line1 + amt.line2 + 
        gene.labels +
        geom_point(data=filter(cell.pairs, !is.na(gene.label) & expr.pass), color='green') +
        geom_point(data=filter(cell.pairs, !is.na(gene.label) & !expr.pass), color='red') +
        scale_color_continuous(name='p-val exp') + 
        scale_size_manual(guide="none",values=c('Other'=.5,'Filtered'=opt.expr.size)) + 
        facet_grid(compare.region~target.region) + xlab(paste(target.label,"\n(log normal mean)")) + ylab(paste(compare.label,"\n(log normal mean)"))

