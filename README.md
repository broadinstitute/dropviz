# dropviz
Prototype app for visualization, exploration of mouse brain single cell gene expression

Data Prep
=========
  
- Create exp_sets.txt from exp_sets.txt.template to identify the names and paths of the experiments to load
- global.R is sourced by the prep- files and also used by shiny. It does a lot of consolidation that could also be cached and loaded.
- run prep-tSNE.R to generate tables of XY coordinates
- run prep-expr.R to generate cached, easily loadable expression data. This takes a long time to run because the DGE matrix is first broken into smaller chunks, then converted into a dense matrix, then converted into a table and then the rows per gene are written to separate files in expr/[exp.label]/gene/[name].RDS. (About 7000 genes (files) per experiment)
XY global coordinates are also stored here because the global coordinates are typically missing or downsampled in plots. Only genes in all.genes are cached at this time, but it can easily be expanded to store all genes.
- run prep-metacells.R to assign the N groups (all but the specified (sub)cluster) and then generate transcript sums and means to be used elsewhere to compute expression levels for scatter plots, differential expression and rankings.
- run prep-marker.R to generate combined marker data (fold change, etc.) in marker/ directory. For all (sub)cluster vs region pairs, store all results per gene. Also the comprehensive gene lists and description, and the top genes. 
- run prep-components.R to generate pre-computed PNGs of components. This is a hack for static images only.

App Logic
=========

The terminology could be improved. The filtering of clusters/subclusters by region/class/cluster/subcluster is referred to as "selected".
Then, for analysis the user choose one cluster or subcluster and one or more genes or components. These are referred to as "current".

Server logic is divided into multiple files with the UI and the model reactives sharing the same file: 

- cell_types.R : the filtering of types (subclusters) by various criteria. provides subclusters.selected() and clusters.selected().
- display_labels.R : the UI allows for regions, subclusters and clusters to be display and queried in different formats (numbers, annotated, all).
  this provides region.names(), cluster.names(), and subcluster.names(), which reactively adjust their .disp column to be used for output
- markers.R : provides cluster.markers() and subcluster.markers() - tables of genes that are differentially expressed according to the current criteria
- tSNE.R : tSNE plots. There is one uber-function, tsne.label, that returns a closure for plotting under the specified parameters
- user_cluster_selection.R : given the selected clusters/subclusters, provides current.cluster() and current.subcluster()
- components.R : independent component - the ICs associated with the current cluster and reactives for returning the cell rotations and XY coordinates.

Cached Data
===========

All pre-computed data is stored in cache/. See the individual prep- files for destination details.
All on-the-fly cached data is stored in www/.
