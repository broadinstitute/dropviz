# dropviz
Shiny app for visualization, exploration of mouse brain single cell gene expression

Data Prep
=========
  
- Create exp_sets.txt from exp_sets.txt.template to identify the names and paths of the experiments to load
- Run prep-all.R
  - prep-tSNE.R to generate tables of XY coordinates
  - prep-expr.R to generate cached, easily loadable expression data. This takes a long time to run because the DGE matrix is first broken into smaller chunks, then converted into a dense matrix, then converted into a table and then the rows per gene are written to separate files in expr/[exp.label]/gene/[name].RDS. (About 7000 genes (files) per experiment)
XY global coordinates are also stored here because the global coordinates are typically missing or downsampled in plots. Only genes in all.genes are cached at this time, but it can easily be expanded to store all genes.
  - prep-metacells.R to assign the N groups (all but the specified (sub)cluster) and then generate transcript sums and means to be used elsewhere to compute expression levels for scatter plots, differential expression and rankings.
  - prep-marker.R to generate combined marker data (fold change, etc.) in marker/ directory. For all (sub)cluster vs region pairs, store all results per gene. Also the comprehensive gene lists and description, and the top genes. 
  - prep-components.R to generate pre-computed PNGs of components. This is a hack for static images only.

Data Files
==========

 <dropviz.prep.dir>/
    expr/
	    <exp>/gene/<gene>.RDS - transcript count, tSNE position and cluster assignment
	metacells/
		<exp>.gene.(sub)clusters.RDS - key <gene, exp, (sub)cluster> expression in (sub)cluster and rest of region, pvalue, fold change, confidence intervals
			intermediate file to create markers/genes/(sub)cluster/<gene>.diffexp.RDS  [prep-marker.R]
    markers/
		genes/(sub)cluster/<gene>.diffexp.RDS - key <exp, (sub)cluster> expression of gene in all (sub)clusters [prep-marker.R]
    tsne
	    [global|local].(sub)clusters.bags.Rdata - bag data
		<exp>/clusterN.xy.RDS - XYs of cells of subclusters in cluster N of exp
		<exp>/global.xy.RDS - XYs of cells in all clusters in exp
		


Cached Data
===========

All pre-computed data is stored in <dropviz.prep.dir>. See the individual prep- files for destination details.
All on-the-fly cached data is stored in www/cache.
