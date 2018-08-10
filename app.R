source("global.R")

server <- function(input, output, session) {

  # display a message, such as a development or testing notice
  if (file.exists("message.txt") && !is.null(getDefaultReactiveDomain())) {
    msg <- readLines("message.txt")
    showNotification(div(h4(msg[1]),msg[2]),duration=NULL,type="warning")
  }

  # majority of the reactive functions.
  source("app-funcs.R", local=TRUE)
  
  # proxy objects for shiny
  source("proxy.R", local=TRUE)

  # output options - do not source when interactive debug
  source("outputOptions.R", local=TRUE)
  
  #####################################################################################################
  # save the latest input for interactive debug
  # load with input <- readRDS("dump.RDS")
  # if bookmarking is enabled, then another option is to simply load the bookmark data, e.g.
  # input <- readRDS("shiny_bookmarks/8628603e8215e159/input.rds")
  observeEvent(input$dump, {
    # deprecated. Use bookmarks
    write.log("Dump"); saveRDS(reactiveValuesToList(input), file="www/dump.RDS")
  })
  
  observeEvent(input$clear.cache, {
    unlink("www/cache/*.png")
  })
  
  # click events on landing page
  observeEvent(input$select.analysis.tab, {
    updateNavbarPage(session, "top-nav", selected = "Query")
  })  
  observeEvent(input$select.go.1, {
    updateNavbarPage(session, "top-nav", selected = "Query")
    updateTabsetPanel(session, "clusterpanel", selected = "tsne")
  })  
  observeEvent(input$select.go.2, {
    updateNavbarPage(session, "top-nav", selected = "Query")
  })  
  observeEvent(input$select.go.3, {
    updateNavbarPage(session, "top-nav", selected = "Query")
  })  
  
  # form reset
  observeEvent(input$resetDisplay, { reset("display-panel") }, ignoreInit = TRUE)
  
  # some special parameters are server-side generated through uiOutput. Calling reset on these has
  # unpredictable, usually undesirable, behavior. So instead, reset a div that doesnot include the
  # uiOutput elements and individually reset other UI elements where necessary.
  observeEvent(input$resetClusters, { 
    updateSelectInput(session, "current.cluster", selected="")
    updateSelectInput(session, "current.subcluster", selected="")  
    reset("compare.multiple")
    reset("cluster-panel") 
  }, ignoreInit = TRUE)
  
  observeEvent(input$resetQuery, {
    filter.vals.init()
    updateSelectInput(session, "cell.type", selected="")
    updateSelectInput(session, "cell.cluster", selected="")
    updateSelectInput(session, "cell.class", selected="")
    updateSelectInput(session, "tissue", selected="")
    reset("showSubclustersInGlobal")
    reset("user.genes")
    updateQueryString("?#")
  }, ignoreInit = TRUE)

  # change the location bar URL when bookmarking and display a pop-up
  onBookmarked(function(url) {
    updateQueryString(url)
    showBookmarkUrlModal(url)
  })

  session$onSessionEnded(trim.cache)
  
  transfer.ics()

}

trim.cache <- function() {
  write.log("Clean cache to 2G of most recently used")
  files <- list.files(cache.dir, full.names=TRUE)

  cache.files <- as_tibble(file.info(files)) %>%
    mutate(files=files) %>%
    filter(!isdir) %>%
    arrange(desc(atime)) %>%
    mutate(csize=cumsum(size))

  old.files <- filter(cache.files, csize>2e9)

  if (nrow(old.files) > 0) {
    print(dplyr::select(cache.files, files, atime))
    unlink(old.files$files)
  }
}

# transfer the IC images, stored in prep.dir, to the cache.dir on first run 
transfer.ics <- function() {
  if (!file.exists(file.path(cache.dir,'ic'))) {
    file.copy(file.path(prep.dir,'ic'), cache.dir, recursive = TRUE)
  }
}

help.doc <- list(tsne.local.label.dl=withTags(span(h4("Help for t-SNE plot of subclusters in local cluster space."),
                                                   p("Each point represents a single cell. Each cell is associated with a gene expression vector. This high-dimensional data within a cluster is reduced using a set of curated independent components and projected onto two dimensions using t-SNE ('local cluster space'). The subcluster classifications are derived from Louvain clustering using a subset of the ICs."),
                                                   p('Subcluster regions are highlighted in different colors based on the filtering choices in the',b('Query'),'section in the left panel. All points in the corresponding cluster are displayed with points outside of the chosen subcluster(s) shown in gray and all points in the chosen subclusters displayed in color. (There is no subsampling for subcluster displays.)'),
                                                   p('If a row in the differentially expressed genes table, below, or one or more genes are entered by name in the',b('Query'), 'panel to the left, then the selected genes\' expression levels will be displayed in black (one row of t-SNE plots per gene) and the mean expression level for that gene in the subcluster is displayed by a color gradient or transparency, depending on settings.'),
                                                   p('Labels and other display features can be customized in the ',b('Display'),'panel'))), 
                 tsne.global.subcluster.label.dl=withTags(span(h4("Help for t-SNE plot of subclusters in global region space."),
                                                               p("This plot is like the t-SNE plot of clusters in global region space, but displays the subcluster labels for each cell.")
                                                               )),  
                 tsne.global.cluster.label.dl=withTags(span(h4("Help for t-SNE plot of clusters in global region space."),
                                                            p("Each point represents a single cell. Each cell is associated with a gene expression vector. This high-dimensional data is reduced using a set of automated independent components and projected onto two dimensions using t-SNE ('global space', i.e. representing all cells from a brain region/tissue). The cluster classifications are derived from Louvain clustering of the ICs."),
                                                            p("Points are generally sub-sampled to improve display and speed rendering. Sampling can be controlled using the ",b('Display'),"panel on the left. 'Bag' plots show the distribution of all points similar to a one-dimensional box plot. The darker region represents 50% of cells. The lighter region represents all points except outliers."),
                                                            p('Clusters are highlighted in different colors based on the filtering choices in the',b('Query'),'section in the left panel. Labels and other display features can be customized in the ',b('Display'),'panel'),
                                                            p('If a row in the differentially expressed genes table, below, or one or more genes are entered by name in the',b('Query'), 'panel to the left, then the selected genes\' expression levels will be displayed in black (one row of t-SNE plots per gene) and the mean expression level for that gene in the subcluster is displayed by a color gradient or transparency, depending on settings.'))), 
                 gene.expr.scatter.subcluster.dl=withTags(span(h4('Subcluster scatter plot'),
                                                               p('Each point is the mean log normalized transcript count among all cells in the target and comparison subclusters (or region). Large points meet the fold ratio and transcript amount criteria in the ',b('Clusters'),'panel. Points are shaded according to their significance. Selected rows in the table of differentially expressed genes are displayed in green or red depending on whether they pass the criteria.'))),
                 gene.expr.scatter.cluster.dl=withTags(span(h4('Cluster scatter plot'),
                                                            p('Each point is the mean log normalized transcript count among all cells in the target and comparison clusters (or region). Large points meet the fold ratio and transcript amount criteria in the ',b('Clusters'),'panel. Points are shaded according to their significance. Selected rows in the table of differentially expressed genes are displayed in green or red depending on whether they pass the criteria.'))),
                 dt.cluster.markers.dl=withTags(span(h4("Differentially expressed genes in clusters"),
                                                     p("Select a 'target' cluster in the left ",b('Clusters')," panel to display those genes that are over-expressed in that cluster with respect to the remaining cells in the region or a chosen comparison cluster. Adjust filter criteria using the ",b("Clusters"),"panel."),
                                                     p("Any manually added genes are always displayed in the table and colored green if the expression criteria is met and colored red otherwise."),
                                                     p("One or more rows can be selected to display gene expression in the t-SNE and scatter plots."))),
                 dt.subcluster.markers.dl=withTags(span(h4("Differentially expressed genes in subclusters"),
                                                        p("Select a 'target' subcluster in the left ",b('Clusters'),"panel to display those genes that are over-expressed in that subcluster with respect to the remaining cells in the region or a chosen comparison subcluster. Adjust filter criteria using the ",b("Clusters"),"panel."),
                                                        p("Any manually added genes are always displayed in the table and colored green if the expression criteria is met and colored red otherwise."),
                                                        p("One or more rows can be selected to display gene expression in the t-SNE and scatter plots, above."))),
                 config=withTags(span(h4("Configuration Panels"),
                                      h5('Query'),p('The configuration panel allows for filtering, comparisons and display adjustments. The ',b('Query'),'panel accepts entry of one or more gene symbols. If entered, then the plots on the right will display the expression with respect to the entered genes. The panel provides auto-complete for a large set of commonly used genes, but other named genes can also be entered.'),
                                      p('In addition to gene entry, the ',b('Query'),'panel allows filtering of the dataset to focus on a subset of the regions, classes, clusters or subclusters. For example, you can limit your search to only "Hippocampus" or to only "Neurons" or to a specific named cluster or subcluster. Choose the t-SNE display on the right to easily visualize what data matches your filtering criteria.'),
                                      h5('Clusters'),p('The ',b('Clusters'),'panel accepts entry of "target" and "comparison" clusters or subclusters (depending on the current plot display). If the meta-group option is enabled, then more than one cluster or subcluster may be combined as the target or comparison. A meta-group sums the expression levels across the chosen clusters or subclusters. Once target and comparisons are selected, then a table is displayed in the bottom right of those genes that are differentially expressed between the two groups. Multiple parameters are provided to modify the classification criteria. Choose the ',b('Scatter'),' panel on the right to visualize how parameter changes affect classification.'),
                                      h5('Display'),p('The ',b('Display'),'panel provides parameters for changing some of the plotting parameters. Only parameters that are relevant to the currently displayed plot on the right are available to the user.'))),
                 gene.expr.rank.cluster.dl=withTags(span(h4("Cluster Levels"),
                                                         p("The plot displays the relative expression of the selected gene(s) among clusters. The error bars represent binomially distributed sampling noise given the number of cells in each cluster; it does not represent heterogeneity in expression among cells within the cluster. Only clusters filtered in the ",b('Query'),"panel are displayed. The plot is limited to only two entered genes. If more than two genes are chosen, then the display switches to a heatmap-like table."))),
                 gene.expr.rank.subcluster.dl=withTags(span(h4("Subcluster Levels"),
                                                            p("The plot displays the relative expression of the selected gene(s) among subclusters. The error bars represent binomially distributed sampling noise given the number of cells in each subcluster; it does not represent heterogeneity in expression among cells within the subcluster. Only subclusters filtered in the ",b('Query'),"panel are displayed. The plot is limited to only two entered genes. If more than two genes are chosen, then the display switches to a heatmap-like table."))),
                 gene.expr.heatmap.cluster.dl=withTags(span(h4("Cluster Levels (Heatmap Table)"),
                                                            p("The table displays the relative expression of the selected gene(s) in all clusters using a shaded dot. Darker dots represent higher gene expression. Only clusters filtered in the ',b('Query'),'panel are displayed. Mouseover dots to display numeric values or choose the download icon to retrieve the full table data in R. The table is displayed when more than two genes are chosen."))),
                 gene.expr.heatmap.subcluster.dl=withTags(span(h4("Subluster Levels (Heatmap Table)"),
                                                               p("The table displays the relative expression of the selected gene(s) in all subclusters using a shaded dot. Darker dots represent higher gene expression. Only subclusters filtered in the ",b('Query'),"panel are displayed. Mouseover dots to display numeric values or choose the download icon to retrieve the full table data in R. The table is displayed when more than two genes are chosen."))),
                 dt.clusters.dl=withTags(span(h4("Filtered Table of Clusters"),
                                              p("The table displays all matching clusters based on the filtering parameters set in the ',b('Query') panel on the left. Other plots and display outputs are limited to the clusters listed here."))),
                 dt.subclusters.dl=withTags(span(h4("Filtered Table of Subclusters"),
                                                 p("The table displays all matching subclusters based on the filtering parameters set in the ',b('Query') panel on the left. Other plots and display outputs are limited to the subclusters listed here.")))
                 )

# draw download icon
downloadIcon <- function(label, title, ...) {
  div(class="top-right",
      downloadLink(label, span(class="glyphicon glyphicon-download-alt" , style="color:#31419a"), 'data-toggle'="tooltip", title=title), ...)
}

helpIcon <- function(label, help.text=NULL) {
  if (is.null(help.text)) {
    help.text <- help.doc[[label]]
  }
  if (!is.null(help.text)) {
    tooltip.span <- span(class="xtooltiptext", help.text) 
  } else {
    tooltip.span <- ''
  }
  
  div(class="xtooltip",
      span(class="glyphicon glyphicon-question-sign", style="color:#31419a"),
      tooltip.span)
}

# draw download and tooltip help icons
downloadAndHelp <- function(label, download.tip="", help.text=NULL) {
  
  downloadIcon(label, title=download.tip, 
               helpIcon(label, help.text))
}

plotDownload <- function(label) {
  downloadAndHelp(label, "Download R Plot")
}
tableDownload <- function(label) {
  downloadAndHelp(label, "Download CSV Table")
}

# replaces <1> in markup with tag[1], <2> with tag[2], etc.
embed.tags <- function(markup, tag) {
  lapply(1:length(tag), function(i) {
    tag.regex <- paste0('<',i,'>')
    markup <<- sub(tag.regex, tag[[i]], markup)  
  })
  markup
}

debug.controls <- function() {
  if (getOption("dropviz.debug", default=FALSE)) {
    div(class="control-box",
        h4("Debug"),
        actionButton("dump","Save State"),
        fluidRow(column(7,checkboxInput("opt.use.cache", "Use Cached Plot Images", value=TRUE)), 
                 column(2,actionButton("clear.cache","Clear Cache"))),
        p(glue("Hostname: {system2('hostname', stdout=TRUE)}")))
  } else {
    span()
  }
}

ui <-
function(request) {
  require(shinyjs)
  
  fluidPage(
    useShinyjs(),
    tags$head(HTML('<script async src="https://www.googletagmanager.com/gtag/js?id=UA-111320548-1"></script>
')),
    tags$head(includeScript("gtag.js")),
    includeCSS("styles.css"),
    tags$link(type="text/css", rel="stylesheet", href="http://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css"),
    navbarPage("DropViz", id="top-nav",
               tabPanel("Home",div(column(2), column(8, embed.tags(HTML(readLines("html/landing.html")), 
                                                                   list(actionButton("select.analysis.tab","Get Started", class="btn btn-lg btn-primary")))), column(2)),
                        embed.tags(HTML(readLines("html/featurette.html")),
                                   list(actionButton("select.go.1", HTML("Explore Cell Clusters &raquo;"), class="btn btn-default"),
                                        actionButton("select.go.2", HTML("Discover Genes &raquo;"), class="btn btn-default"),
                                        actionButton("select.go.3", HTML("Search &raquo;"), class="btn btn-default")))),
               tabPanel("Query",
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(width=3, id="controlpanel",
                                       h2('Parameters', class="Param"),
                                       div(helpIcon("config"), class='top-right'),
                                       tabsetPanel(type="tabs", id='controltabs', selected="controltabs-query",
                                                   tabPanel("Query", value="controltabs-query",
                                                            div(id="filter-params", class="control-box", style="margin-bottom: 0",
                                                                fluidRow(
                                                                  column(12,
                                                                         selectizeInput("user.genes", "Gene", choices=c("Symbol"=""),
                                                                                        multiple=TRUE, width='100%', options=list(create=TRUE,persist=FALSE,closeAfterSelect=TRUE))
                                                                         )
                                                                ),
                                                                fluidRow(
                                                                  column(12, uiOutput("region"))),
                                                                fluidRow(column(12, uiOutput("cell.class"))),
                                                                fluidRow(column(12, uiOutput("cell.cluster"))), 
                                                                conditionalPanel("input.mainpanel=='subclusters'",
                                                                                 fluidRow(column(12, uiOutput("cell.type")))
                                                                ),
                                                                conditionalPanel("input.mainpanel=='subclusters' && input.subclusterpanel=='tsne'",
                                                                                 checkboxInput("showSubclustersInGlobal","Show Subclusters in Global Plot", value=FALSE)),
                                                                fluidRow(style="margin-bottom: 25px",
                                                                         column(5, offset=1, actionButton("go", "Update!")),
                                                                         column(6, actionButton("resetQuery", "Reset Query")))
                                                            )
                                                   ),
                                                   tabPanel("Clusters", value="controltabs-compare",
                                                            div(class="control-box",
                                                                conditionalPanel('input.mainpanel=="clusters"',
                                                                                 h4("Compare Clusters"),
                                                                                 fluidRow(
                                                                                   uiOutput("current.cluster"),
                                                                                   uiOutput("comparison.cluster")
                                                                                 )
                                                                                 ),
                                                                conditionalPanel('input.mainpanel=="subclusters"',
                                                                                 h4("Compare Subclusters"),
                                                                                 fluidRow(
                                                                                   uiOutput("current.subcluster"),
                                                                                   uiOutput("comparison.subcluster")
                                                                                 )
                                                                                 ),
                                                                checkboxInput("compare.multiple","Allow target/comparison meta-groups", value=FALSE)
                                                                ),
                                                            div(id="cluster-panel",
                                                                div(class="control-box",
                                                                    h4("Differential Expression Criteria"),
                                                                    sliderInput("fold.change", "Minimum Fold Ratio", min=1, max=30, value=2, step=0.5),
                                                                    sliderInput("pval.thresh", "Maximum P-Value Exponent", min=-300, max=0, value=-100, step=1),
                                                                    selectInput("expr.filter.opt", "", c("AND"="both","OR"="either","Fold Ratio ONLY (two-sided)"="fc","Abundance ONLY"="amt")),
                                                                    tags$table(tags$tr(tags$td(sliderInput("min.amt.within", "Min Mean Log Amount in Target", min=0, max=6, value=2.5, step=0.25),
                                                                                               valign="top"),
                                                                                       tags$td(width="5%"),
                                                                                       tags$td(sliderInput("max.amt.without", "Max Mean Log Amount in Comp", min=0, max=6, value=1, step=0.25),
                                                                                               valign="top")),
                                                                               width='100%'),
                                                                    div(style="display:none", actionButton("upload.genes","Upload Gene List", width='100%', onclick="alert('Not Implemented')"))
                                                                ),
                                                                fluidRow(style="margin-top: 10px",
                                                                         column(10, offset=2, actionButton("resetClusters", "Reset Cluster Parameters")))
                                                            )),
                                                   tabPanel("Display",
                                                            div(id='display-panel',
                                                                conditionalPanel('(input["mainpanel"]=="clusters" && input["clusterpanel"]=="rank") || (input["mainpanel"]=="subclusters" && input["subclusterpanel"]=="rank")',
                                                                                 conditionalPanel('input["user.genes"]==undefined || input["user.genes"].length <= 2',
                                                                                                  div(class="control-box", h4("Level Plot Settings"), 
                                                                                                      checkboxInput("opt.rank.by.region", "Group Rankings By Region", value=FALSE), 
                                                                                                      selectInput("top.N","Top N Cluster or Subcluster", choices=c(1,2,3,4,5,10,20,50,200),selected=10),
                                                                                                      selectInput("opt.plot.height", "Plot Height", choices=c('Fixed'='fixed', 'Adjust'='auto'), selected='auto')
                                                                                                  )),
                                                                                 conditionalPanel('input["user.genes"]!=undefined && input["user.genes"].length > 2',
                                                                                                  div(class="control-box", h4("Heat Map Settings"), 
                                                                                                      sliderInput("opt.heatmap.max.per100k", "Threshold Max Transcripts (per 100k) in Heatmap Display (lower values increase sensitivity)",
                                                                                                                  0, 1000, value=100, step=10)
                                                                                                  ))
                                                                ),
                                                                conditionalPanel("input['user.genes'] && ((input.mainpanel=='clusters' && input.clusterpanel=='tsne') || (input.mainpanel=='subclusters' && input.subclusterpanel=='tsne'))",
                                                                                 div(class="control-box",
                                                                                     h4("Gene Search Settings"),
                                                                                     div(style="display:none",checkboxInput("normalize.expression.by.facet","Normalize Expression within Region or Cluster", value=FALSE)),
                                                                                     selectizeInput("opt.tx", "Show Metacell Expression as", choices=c("Overlay"="alpha", "Hot-Cool"="heat"), selected = "alpha"),
                                                                                     sliderInput("opt.tx.min", "Show Labels When Metacell Expression is Greater than % of Max", 0, 90, value=70, round=TRUE, step=10, post='%'),
                                                                                     conditionalPanel("input['opt.tx']=='heat'",
                                                                                                      checkboxInput("opt.tx.legend","Show Metacell Expression Legend", value=FALSE)),
                                                                                     selectInput("opt.tx.scale", "Scale transparency or color range: ", choices=c("Fixed"="fixed", "Observed Max Value"="max"), selected="max"),
                                                                                     checkboxInput("opt.tx.cells","Show Expression Per Cell", value=TRUE),
                                                                                     div(style="display:none",conditionalPanel("!input['opt.tx.cells'] && input['user.genes'].length > 1", checkboxInput("opt.tx.sum", "Sum Expression of Multiple Search Genes", value=FALSE))))
                                                                ),
                                                                conditionalPanel('(input["mainpanel"]=="clusters" && input["clusterpanel"]=="tsne") || (input["mainpanel"]=="subclusters" && input["subclusterpanel"]=="tsne")',
                                                                                 conditionalPanel('input["opt.tx.cells"]',
                                                                                                  div(class="control-box",
                                                                                                      h4("Cell Expression Settings"),
                                                                                                      selectInput("opt.cell.display.type","Display Gene Expression Using", choices=c("Size"="size","Absent/Present"="detect"), selected="size"),
                                                                                                      conditionalPanel('input["opt.cell.display.type"]=="detect"',
                                                                                                                       div(style="display:none",sliderInput("opt.detection.thresh","Observed Transcript Copies (FIXME - currently data is normalized values)", 0, 5, value=0, step=1))),
                                                                                                      conditionalPanel("input['opt.cell.display.type']!='detect'",
                                                                                                                       sliderInput("opt.expr.size", "Point Size of Maximum Expression", 1, 10, value=4, step=1))
                                                                                                  )
                                                                                 ),
                                                                                 div(class="control-box",
                                                                                     h4("t-SNE Plot Settings"),
                                                                                     selectInput("opt.plot.label", "Plot Labels for Clusters and Subclusters", choices=c("Names"='disp',"Numbers"='number',"None"='none')),
                                                                                     conditionalPanel("!input['user.genes'] || (input['opt.tx']=='alpha' && ((input.mainpanel=='clusters' && input.clusterpanel=='tsne') || (input.mainpanel=='subclusters' && input.subclusterpanel=='tsne')))",
                                                                                                      checkboxInput("opt.show.bags","Display t-SNE using bag plots", value=FALSE)),
                                                                                     selectInput("opt.downsampling.method","Downsample Cells",choices=c("Uniformly"='uniform',"Per Cluster"='cluster',"Show all"='none'), selected='uniform'),
                                                                                     conditionalPanel("input['opt.downsampling.method']!='none'",
                                                                                                      sliderInput("downsampling", "Downsample Count", 0, 100000, value=6000, step=1000))
                                                                                 )
                                                                ),
                                                                conditionalPanel('(input["mainpanel"]=="clusters" && input["clusterpanel"]=="scatter") || (input["mainpanel"]=="subclusters" && input["subclusterpanel"]=="scatter")',
                                                                                 div(class="control-box", 
                                                                                     h4("Scatter Plot Settings"),
                                                                                     checkboxInput("opt.scatter.gene.labels","Show Gene Labels on Scatter Plots", value=TRUE)
                                                                                 )
                                                                ),
                                                                div(class="control-box",
                                                                    h4("Labels"),
                                                                    div(style="display:none", selectInput("opt.cluster.disp","Label Clusters", choices=c("Using Annotated Class and Markers"='annotated', "With Numbers"='numbers',"Class, Markers and Numbers"="all"), selected = 'all')),
                                                                    div(style="display:none", selectInput("opt.region.disp","Label Region", choices=c("Using Region Name"='region',"Using Experiment Name"='experiment'))),
                                                                    conditionalPanel("input['opt.cluster.disp']=='annotated' || input['opt.cluster.disp']=='all'",
                                                                                     checkboxInput("use.common.name", "Use Common Name for Subcluster, If Present", value = FALSE)),
                                                                    conditionalPanel("input['use.common.name']", span(h6("(Common names are interpretive best guesses)")))),
                                                                fluidRow(style="margin-top: 10px", column(10, offset=2, actionButton("resetDisplay", "Reset Display Parameters")))
                                                                )
                                       ),
                                       fluidRow(style="margin-top:20px", column(10, offset=2, bookmarkButton("Bookmark / Update URL"))),
                                       debug.controls()
                          )),
                          
                          # Show a plot of the generated distribution
                          mainPanel(width=9, 
                                    tabsetPanel(type="tabs", id="mainpanel",
                                                tabPanel("Global Clusters", value = "clusters",
                                                         tabsetPanel(type="pills", id="clusterpanel",
                                                                     tabPanel("Levels By Cluster", value="rank",
                                                                              div(conditionalPanel("input['user.genes']==undefined || input['user.genes'].length <= 2",
                                                                                                   plotDownload("gene.expr.rank.cluster.dl"),
                                                                                                   uiOutput("gene.expr.rank.cluster.output")
                                                                                                   ),
                                                                                  conditionalPanel("input['user.genes']!=undefined && input['user.genes'].length > 2",
                                                                                                   tableDownload("gene.expr.heatmap.cluster.dl"),
                                                                                                   DT::dataTableOutput("gene.expr.heatmap.cluster")
                                                                                                   ))),
                                                                     tabPanel("tSNE", value="tsne",
                                                                              plotDownload("tsne.global.cluster.label.dl"),
                                                                              fluidRow(div(id="global-expression", class="scroll-area", 
                                                                                           span(class="img-center",imageOutput("tsne.global.cluster.label",height="500px"))))),
                                                                     tabPanel("Scatter", value="scatter",
                                                                              plotDownload("gene.expr.scatter.cluster.dl"),
                                                                              fluidRow(div(id="global-scatter", class="scroll-area",
                                                                                           span(class="img-center",imageOutput("gene.expr.scatter.cluster", height=500))))),
                                                                     tabPanel("Table",
                                                                              tableDownload("dt.clusters.dl"),
                                                                              fluidRow(div(id="dt-clusters", class="scroll-area", 
                                                                                           DT::dataTableOutput("dt.clusters"))))),
                                                         hr(),
                                                         uiOutput("dt.cluster.markers.heading"),
                                                         conditionalPanel(
                                                           'input["current.cluster"]',
                                                           fluidRow(class="table-area",
                                                                    tableDownload("dt.cluster.markers.dl"),
                                                                    column(11,DT::dataTableOutput("dt.cluster.markers")))
                                                         )
                                                ),
                                                tabPanel("Subclusters", value = "subclusters",
                                                         tabsetPanel(type="pills", id="subclusterpanel",
                                                                     tabPanel("Levels By Subcluster", value="rank",
                                                                              div(conditionalPanel("input['user.genes']==undefined || input['user.genes'].length <= 2",
                                                                                                   plotDownload("gene.expr.rank.subcluster.dl"),
                                                                                                   uiOutput("gene.expr.rank.subcluster.output")
                                                                                                   ),
                                                                                  conditionalPanel("input['user.genes']!=undefined && input['user.genes'].length > 2",
                                                                                                   tableDownload("gene.expr.heatmap.subcluster.dl"),
                                                                                                   DT::dataTableOutput("gene.expr.heatmap.subcluster")
                                                                                                   ))),
                                                                     tabPanel("tSNE", value="tsne",
                                                                              conditionalPanel("input.showSubclustersInGlobal",
                                                                                               plotDownload("tsne.global.subcluster.label.dl")),
                                                                              conditionalPanel("!input.showSubclustersInGlobal",
                                                                                               plotDownload("tsne.local.label.dl")),
                                                                              fluidRow(div(id="local-expression", class="scroll-area",
                                                                                           conditionalPanel("input.showSubclustersInGlobal",
                                                                                                            span(class="img-center",imageOutput("tsne.global.subcluster.label", height=500))),
                                                                                           conditionalPanel("!input.showSubclustersInGlobal",
                                                                                                            span(class="img-center",imageOutput("tsne.local.label", height=500)))))),
                                                                     tabPanel("Scatter", value="scatter",
                                                                              plotDownload("gene.expr.scatter.subcluster.dl"),
                                                                              fluidRow(div(id="local-scatter", class="scroll-area",
                                                                                           span(class="img-center",imageOutput("gene.expr.scatter.subcluster", height=500))))),
                                                                     tabPanel("Table",
                                                                              tableDownload("dt.subclusters.dl"),
                                                                              fluidRow(div(id="dt-subclusters", class="scroll-area", DT::dataTableOutput("dt.subclusters"))))),
                                                         hr(),
                                                         conditionalPanel('input["subclusterpanel"]!="IC"',
                                                                          uiOutput("dt.subcluster.markers.heading"),
                                                                          conditionalPanel(
                                                                            'input["current.subcluster"]',
                                                                            fluidRow(class="table-area",
                                                                                     tableDownload("dt.subcluster.markers.dl"),
                                                                                     column(11,DT::dataTableOutput("dt.subcluster.markers")))))
                                                         ),
                                                tabPanel("Independent Components", value="IC",
                                                         fluidRow(div(id="ic-grid", class="scroll-area",
                                                                      uiOutput("ic.grid"))),
                                                         hr(),
                                                         uiOutput("dt.components.heading"),
                                                         fluidRow(DT::dataTableOutput("dt.components"))
                                                         )
                                                )
                                    )
                        )
               ),
               tabPanel("Team",
                        HTML(readLines("html/team.html"))),
               tabPanel("Data",
                        h2("Data Downloads"),hr(),
                        tags$dl(tags$dt("Metacells"),
                                tags$dd(p("Gene expression profiles of the 565 transcriptionally distinct cell populations identified across nine regions in the adult mouse brain. Each column is a \"metacell.\" There is one metacell for every subcluster, which contains the aggregate UMI counts for all the single-cells that belong to that subcluster."),
                                        p(a(href="https://storage.googleapis.com/dropviz-downloads/static/metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv","CSV: metacells.BrainCellAtlas_Saunders_version_2018.04.01.csv (40M)")),
                                        p(a(href="https://storage.googleapis.com/dropviz-downloads/static/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS","R Data: metacells.BrainCellAtlas_Saunders_version_2018.04.01.rds (15M)"))),
                                tags$dt("Annotations"),
                                tags$dd(p("Annotation file for the 565 atlas cell populations. Provides the tissue of origin, cell class, formal markers, formal full name and common name (anatomical best guess) for each metacell."),
                                        p(a(href="https://storage.googleapis.com/dropviz-downloads/static/annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv","CSV: annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv (54K)")),
                                        p(a(href="https://storage.googleapis.com/dropviz-downloads/static/annotation.BrainCellAtlas_Saunders_version_2018.04.01.xlsx","Excel: annotation.BrainCellAtlas_Saunders_version_2018.04.01.xlsx (83K)")),
                                        p(a(href="https://storage.googleapis.com/dropviz-downloads/static/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS","R Data: annotation.BrainCellAtlas_Saunders_version_2018.04.01.rds (12K)"))),
                                tags$dt("Single Cell Suspension Protocol from Acute Adult Brain"),
                                tags$dd(p(),p(a(href="https://storage.googleapis.com/dropviz-downloads/static/Saunders_scBrainSuspensionProtocol_v1_180419.pdf","Saunders_scBrainSuspensionProtocol_v1_180419.pdf"))),
                                tags$dt("Instruction for loading DGE files"),
                                HTML(readLines("html/dge-instructions.html")),
                                tags$dt("DGE By Region"),
                                HTML(readLines("html/dge-by-region.html")),
                                tags$dt("DGE By Class"),
                                HTML(readLines("html/dge-by-class.html")))
                        ),
               tabPanel("Tutorial",
                        h1("Video Tutorials"), hr(),
                        div(style="text-align: center;width: 100%;vertical-align: middle;height: 100%;",
                            tags$iframe(width="640", height="385", src="https://www.youtube.com/embed/videoseries?list=PLRL1w4IELl9j0SuViEj3Zj6HmsBULZbbs",
                                        frameborder="0", allow="autoplay; encrypted-media", allowfullscreen="true"))
                        ),
               tabPanel("Feedback",
                        h3("Feedback"),
                        p("We welcome any comments, bug reports, and feature requests. Please send all feedback to ",a(href="mailto:mouse.dropviz@gmail.com","mouse.dropviz@gmail.com"),"."))
    ),
    tags$script(src="http://code.jquery.com/ui/1.12.1/jquery-ui.min.js")## ,
    ## tags$script("jQuery(function (){ $('.scroll-area').resizable(); });")
  )
}


shinyApp(ui, server, enableBookmarking = "server")

