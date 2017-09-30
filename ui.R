
# this is a hack to set the ID for the actual sidebar div because shiny::sidebarPanel() assigns the
# id to the sidebar's child (a form).
# jsCode <- "
# shinyjs.init = function() { $('#sidebarform').parent().attr('id', 'sidebar') }
# "

help.doc <- list(tsne.local.label.dl=withTags(span(h4("Help for t-SNE plot of subclusters in local cluster space."),
                                                   p("Each point represents a single cell. All points are displayed. Each cell is associated with a gene expression vector. This high-dimensional data within a cluster is reduced using a set of curated independent components and projected onto two dimensions using t-SNE ('local cluster space'). The subcluster classifications are derived from Louvain clustering using a subset of the ICs. (See paper...)"),
                                                   p('Subcluster regions are highlighted in different colors based on the filtering choices in the',b('Highlight'),'section in the left panel.'),
                                                   p('If differentially expressed genes are displayed in the bottom table, then selecting rows in the table will overlay the selected genes\' expression levels with larger size points representing higher transcript count.'),
                                                   p('If independent components are displayed in the bottom table, then selecting rows will overlay the selected ICs with weights displayed in a gradient color (blue: negative, red: positive).'),
                                                   p('If ICs are selected, then any colored subcluster selection is disabled because it is not possible to display both a categorical color (subclusters) and a gradient (IC weight) on the same plot.'))), 
                 tsne.global.subcluster.label.dl=withTags(span(h4("Help for t-SNE plot of subclusters in global region space."),
                                                               p("This plot is similar to the t-SNE plot of clusters in global region space, but displays the subcluster labels for each cell.")
                 )),  
                 tsne.global.cluster.label.dl=withTags(span(h4("Help for t-SNE plot of clusters in global region space."),
                                                            p("Each point represents a single cell. Each cell is associated with a gene expression vector. This high-dimensional data is reduced using a set of automated independent components and projected onto two dimensions using t-SNE ('global space', i.e. representing all cells from a brain region). The cluster classifications are derived from Louvain clustering of the ICs. (See paper...)"),
                                                            p("Points are generally sub-sampled to improve display and speed rendering. Sampling can be controlled using the ",b('Display'),"panel on the left. 'Bag' plots show the distribution of all points (no sampling). The darker region represents 50% of cells. The lighter region represents all points except outliers."),
                                                            p('Clusters are highlighted in different colors based on the filtering choices in the',b('Highlight'),'section in the left panel. Text, numeric or no labels can be set in the ',b('Display'),'panel'),
                                                            p('If differentially expressed genes are displayed in the bottom table, then selecting rows in the table will overlay the selected genes\' expression levels with larger size points representing higher transcript count.')
                 )), 
                 gene.expr.scatter.subcluster.dl=withTags(span(h4('Subcluster scatter plot'),
                                                               p('Each point is the mean log normalized transcript count among all cells in the target and comparison subclusters (or region). Large points meet the fold ratio and transcript amount criteria in the ',b('Genes'),'panel. Points are shaded according to their significance. Selected rows in the table are displayed in green or red depending on whether they pass the criteria.'))),
                 gene.expr.scatter.cluster.dl=withTags(span(h4('Cluster scatter plot'),
                                                            p('Each point is the mean log normalized transcript count among all cells in the target and comparison subclusters (or region). Large points meet the fold ratio and transcript amount criteria in the ',b('Genes'),'panel. Points are shaded according to their significance. Selected rows in the table are displayed in green or red depending on whether they pass the criteria.'))),
                 dt.cluster.markers.dl=withTags(span(h4("Differentially over-expressed genes in clusters."),
                                                     p("Select a 'target' cluster in the left ",b('Cells'),"panel to display those genes that are over-expressed in that cluster with respect to the remaining cells in the region or a chosen comparison cluster. Adjust filter criteria using the ",b("Genes"),"panel on the left."),
                                                     p("Any manually added genes are always displayed in the table and colored green if the expression criteria is met and red otherwise."),
                                                     p("One or more rows can be selected to display gene expression in the t-SNE and scatter plots, above."))),
                 dt.subcluster.markers.dl=withTags(span(h4("Differentially over-expressed genes in subclusters."),
                                                        p("Select a 'target' subcluster in the left ",b('Cells'),"panel to display those genes that are over-expressed in that subcluster with respect to the remaining cells in the region or a chosen comparison subcluster. Adjust filter criteria using the ",b("Genes"),"panel on the left."),
                                                        p("Any manually added genes are always displayed in the table and colored green if the expression criteria is met and red otherwise."),
                                                        p("One or more rows can be selected to display gene expression in the t-SNE and scatter plots, above."))),
                 config='TBD. Help for the left configuration panel'
                 )

# draw download icon
downloadIcon <- function(label, title, ...) {
  div(class="top-right",
      downloadLink(label, span(class="glyphicon glyphicon-download-alt"), 'data-toggle'="tooltip", title=title), ...)
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
      span(class="glyphicon glyphicon-question-sign", style="color:#428BCA"),
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
    tabPanel("X",
             div(class="control-box",
                 h4("Debug"),
                 actionButton("dump","Save State"),
                 fluidRow(column(7,checkboxInput("opt.use.cache", "Use Cached Plot Images", value=TRUE)), 
                          column(2,actionButton("clear.cache","Clear Cache"))),
                 p(glue("Hostname: {system2('hostname', stdout=TRUE)}")))
    )
  } else {
    span()
  }
}

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    useShinyjs(),
    #    extendShinyjs(text = jsCode),
    includeCSS("styles.css"),
    tags$link(type="text/css", rel="stylesheet", href="http://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css"),
    navbarPage("DropViz", id="top-nav",
               tabPanel("Home",div(column(2), column(8, embed.tags(HTML(readLines("html/landing.html")), 
                                                                   list(actionButton("select.analysis.tab","Get Started", class="btn btn-lg btn-primary")))), column(2)),
                        HTML(readLines("html/featurette.html"))),
               tabPanel("Analysis",
                        
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(width=3, id="controlpanel",
                                       div(helpIcon("config"), class='top-right'),
                                       tabsetPanel(type="tabs",
                                                   tabPanel("Query",
                                                            div(class="control-box", style="margin-bottom: 0",
                                                                h4("Filter"),
                                                                fluidRow(
                                                                  column(6, selectizeInput("user.genes", "Gene", choices=c("Symbol"="",top.genes),
                                                                                           multiple=TRUE, width='100%', options=list(create=TRUE,persist=FALSE))),
                                                                  column(6, uiOutput("region"))),
                                                                conditionalPanel("input['user.genes'] && input['opt.show.bags'] && ((input.mainpanel=='clusters' && input.clusterpanel=='tSNE') || (input.mainpanel=='subclusters' && input.subclusterpanel=='tSNE'))",
                                                                                 checkboxInput("opt.tx.alpha",span(style="font-size: small","Show Expression as Transparency in t-SNE"),value = TRUE)),
                                                                conditionalPanel("input.mainpanel=='subclusters' && input.subclusterpanel=='tSNE'",
                                                                                 checkboxInput("showSubclustersInGlobal","Show Subclusters in Global Plot", value=FALSE)),
                                                                fluidRow(
                                                                  column(6, uiOutput("cell.class")),
                                                                  column(6, uiOutput("cell.cluster")), 
                                                                  conditionalPanel("input.mainpanel=='subclusters'",
                                                                                   column(12, uiOutput("cell.type")))
                                                                )
                                                            )
                                                   ),
                                                   tabPanel("Compare",
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
                                                                )
                                                            ),
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
                                                                actionButton("upload.genes","Upload Gene List", width='100%', onclick="alert('Not Implemented')")
                                                            )
                                                   ),
                                                   tabPanel("Display",
                                                            # opt.*.label effects how items are displayed in pull-downs, table titles, etc. If these are changed
                                                            # in real time, then options generally reset. One advantage of different displayed options is that querying
                                                            # can be more general. E.g. If cluster labels are just a class and marker e.g. Atrocyte.Gja1 instead of
                                                            # Astrocyte.Gja1 [#4], then you can match the same cluster across different regions. But if these are
                                                            # changed on the fly, then current selections are lost. And there's no general solution because there
                                                            # isn't a 1:1 mapping of display options to model values.
                                                            div(style="display:none",
                                                                selectInput("opt.cluster.disp","Label Clusters", choices=c("Using Annotated Class and Markers"='annotated', "With Numbers"='numbers',"Class, Markers and Numbers"="all"), selected = 'all'),
                                                                selectInput("opt.region.disp","Label Region", choices=c("Using Region Name"='region',"Using Experiment Name"='experiment')),
                                                                conditionalPanel("input['opt.cluster.disp']=='annotated' || input['opt.cluster.disp']=='all'",
                                                                                 checkboxInput("use.common.name", "Use Common Name for Subcluster, If Present", value = TRUE))),
                                                            selectInput("opt.plot.label", "Plot Labels for Clusters and Subclusters", choices=c("Names"='disp',"Numbers"='number',"None"='none')),
                                                            checkboxInput("opt.show.bags","Display t-SNE using bag plots", value=TRUE),
                                                            selectInput("opt.downsampling.method","Downsample Cells",choices=c("Uniformly"='uniform',"Per Cluster"='cluster',"Show all"='none'), selected='uniform'),
                                                            conditionalPanel("input['opt.downsampling.method']!='none'",
                                                                             sliderInput("downsampling", "Downsample Count", 0, 100000, value=2000, step=1000)),
                                                            selectInput("top.N","Gene Search: Show Top Cluster or Subcluster Matches", choices=c(1,2,3,4,5,10,20),selected=5),
                                                   
                                                            sliderInput("opt.expr.size", "Point Size of Maximum Expression", 1, 10, value=3, step=0.5),
                                                            checkboxInput("opt.scatter.gene.labels","Show Gene Labels on Scatter Plots", value=TRUE),
                                                            selectInput("opt.components", "Show Components", choices=c("Real"='real','Used for Clustering'='clustering','All'='all'))
                                                   ),
                                                   debug.controls()
                                       )
                          ),
                          
                          # Show a plot of the generated distribution
                          mainPanel(width=9, 
                                    tabsetPanel(type="tabs", id="mainpanel",
                                                tabPanel(span("Global Clusters"), value = "clusters",
                                                         tabsetPanel(type="pills", id="clusterpanel",
                                                                     tabPanel("tSNE",
                                                                              fluidRow(div(id="global-expression", class="scroll-area", 
                                                                                           plotDownload("tsne.global.cluster.label.dl"),
                                                                                           span(class="img-center",imageOutput("tsne.global.cluster.label", height="500px"))))),
                                                                     tabPanel("Rank", 
                                                                              fluidRow(div(id="global-rank", class="scroll-area",
                                                                                           plotDownload("gene.expr.rank.cluster.dl"),
                                                                                           span(class="img-center",plotOutput("gene.expr.rank.cluster", height=500))))),
                                                                     tabPanel("Scatter", 
                                                                              fluidRow(div(id="global-scatter", class="scroll-area",
                                                                                           plotDownload("gene.expr.scatter.cluster.dl"),
                                                                                           span(class="img-center",imageOutput("gene.expr.scatter.cluster", height=500))))),
                                                                     tabPanel("Table",
                                                                              fluidRow(div(id="dt-clusters", class="scroll-area", 
                                                                                           DT::dataTableOutput("dt.clusters"))))),
                                                         hr(),
                                                         uiOutput("dt.cluster.markers.heading"),
                                                         conditionalPanel(
                                                             'true && !$("div[id=\'current.cluster\']").is(":empty")',
                                                             fluidRow(class="table-area",
                                                                      tableDownload("dt.cluster.markers.dl"),
                                                                      column(11,DT::dataTableOutput("dt.cluster.markers")))
                                                         )
                                                ),
                                                tabPanel("Subclusters", value = "subclusters",
                                                         tabsetPanel(type="pills", id="subclusterpanel",
                                                                     tabPanel("tSNE",
                                                                              fluidRow(div(id="local-expression", class="scroll-area",
                                                                                           conditionalPanel("input.showSubclustersInGlobal",
                                                                                                            plotDownload("tsne.global.subcluster.label.dl"),
                                                                                                            span(class="img-center",imageOutput("tsne.global.subcluster.label", height=500))),
                                                                                           conditionalPanel("!input.showSubclustersInGlobal",
                                                                                                            plotDownload("tsne.local.label.dl"),
                                                                                                            span(class="img-center",imageOutput("tsne.local.label", height=500)))))),
                                                                     tabPanel("Rank", 
                                                                              fluidRow(div(id="local-rank", class="scroll-area",
                                                                                           plotDownload("gene.expr.rank.subcluster.dl"),
                                                                                           imageOutput("gene.expr.rank.subcluster", height=500)))),
                                                                     tabPanel("Scatter",
                                                                              fluidRow(div(id="local-scatter", class="scroll-area",
                                                                                           plotDownload("gene.expr.scatter.subcluster.dl"),
                                                                                           span(class="img-center",imageOutput("gene.expr.scatter.subcluster", height=500))))),
                                                                     tabPanel("Independent Components",
                                                                              fluidRow(div(id="ic-grid", class="scroll-area",
#                                                                                           plotDownload("ic.grid.dl"),
                                                                                           span(class="img-center",plotOutput("ic.grid", height=500))))),
                                                                     tabPanel("Table",  
                                                                              fluidRow(div(id="dt-subclusters", class="scroll-area", DT::dataTableOutput("dt.subclusters"))))),
                                                         hr(),
                                                         tabsetPanel(type="pills",
                                                                     tabPanel("Differential Expression",
                                                                              uiOutput("dt.subcluster.markers.heading"),
                                                                              fluidRow(class="table-area",
                                                                                       tableDownload("dt.subcluster.markers.dl"),
                                                                                       column(11,DT::dataTableOutput("dt.subcluster.markers")))),
                                                                     tabPanel("Independent Components",
                                                                              uiOutput("dt.components.heading"),
                                                                              fluidRow(DT::dataTableOutput("dt.components")))))
                                                # ,
                                                # tabPanel("Metacells")
                                    )
                          )
                          #    tabPanel("Community Annotations", p("A wiki-like editable annotation of cell type"))
                        )
               ),
               tabPanel("Data",
                        p("Data sets will be released concurrently with the publication of our manuscript. Please check back soon.")),
               tabPanel("Credits",
                        h3("Contributors"),
                        p("A list of people and perhaps their roles"),
                        h3("Contact"), p(a(href="mailto:asaunders@genetics.med.harvard.edu","A. Saunders"),"and",p(href="mailto:emacosko@broadinstitute.org","E. Macosko")),
                        h3("Funding"), p("A.S is funded through a post-doctoral fellowship from the Helen Hay Whitney foundation."))
    ),
    tags$script(src="http://code.jquery.com/ui/1.12.1/jquery-ui.min.js"),
    tags$script("jQuery(function (){ $('.scroll-area').resizable(); });")
  )
)
