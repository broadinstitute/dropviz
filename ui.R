
# this is a hack to set the ID for the actual sidebar div because shiny::sidebarPanel() assigns the
# id to the sidebar's child (a form).
# jsCode <- "
# shinyjs.init = function() { $('#sidebarform').parent().attr('id', 'sidebar') }
# "

help.doc <- list(tsne.local.label.dl='Help for tSNE plot of subclusters in local tSNE space', 
                 tsne.global.subcluster.label.dl='Help for tSNE plot of subclusters in global tSNE space', 
                 tsne.global.cluster.label.dl='Help for tSNE plot of clusters in global tSNE space',
                 gene.expr.scatter.subcluster.dl='Subcluster scatter plot help goes here!',
                 gene.expr.scatter.cluster.dl='Help for cluster scatter plot',
                 dt.cluster.markers.dl='This is the help text for differentially expressed genes among clusters',
                 dt.subcluster.markers.dl='This is the help text for differentially expressed genes among SUBclusters',
                 config='Help for the left configuration panel'
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

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    useShinyjs(),
    #    extendShinyjs(text = jsCode),
    includeCSS("styles.css"),
    tags$link(type="text/css", rel="stylesheet", href="http://code.jquery.com/ui/1.12.1/themes/smoothness/jquery-ui.css"),
    navbarPage(span("DropViz",br(),"Single Cell Sequencing"),
               # titlePanel("DropViz Prototype"),
               # h4("Cell Types Defined via Large Scale Single Cell Mouse Brain Gene Expression"),
               tabPanel("Clustering",
                        
                        
                        # Sidebar with a slider input for number of bins 
                        sidebarLayout(
                          sidebarPanel(width=3, id="sidebarform",
                                       div(helpIcon("config"), class='top-right'),
                                       tabsetPanel(type="tabs",
                                                   tabPanel("Cells",
                                                            tags$table(tags$tr(uiOutput("region", container=tags$td, width='50%'),
                                                                               uiOutput("cell.class", container=tags$td, width='50%')),
                                                                       tags$tr(uiOutput("cell.cluster", container=tags$td, colspan="2")),
                                                                       tags$tr(uiOutput("cell.type", container=tags$td, colspan="2")),
                                                                       width="100%")
                                                   ),
                                                   tabPanel("Genes",
                                                            h4("Differential Expression Criteria"),
                                                            sliderInput("fold.change", "Minimum Fold Change", min=0, max=30, value=2, step=0.5),
                                                            selectInput("expr.filter.opt", "", c("AND"="both","OR"="either","Fold Change ONLY"="fc","Percent Present ONLY"="pp")),
                                                            tags$table(tags$tr(tags$td(sliderInput("min.pct.within", "Min Percent Present in Target", min=0, max=1, value=.75, step=0.05),
                                                                                       valign="top"),
                                                                               tags$td(width="5%"),
                                                                               tags$td(sliderInput("max.pct.without", "Max Percent Present in Comparison", min=0, max=1, value=.25, step=0.05),
                                                                                       valign="top")),
                                                                       width='100%'),
                                                            hr(),
                                                            h4("Manual Gene Selection"),
                                                            selectizeInput("user.genes", "", choices=c("Choose Additional Genes Of Interest"="",all.genes),
                                                                           multiple=TRUE, width='100%'),
                                                            actionButton("upload.genes","Upload Gene List", width='100%', onclick="alert('Not Implemented')")
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
                                                            selectInput("opt.downsampling.method","Down Sample Cells",choices=c("Uniformly"='uniform',"Per Cluster"='cluster',"Show all"='none'), selected='uniform'),
                                                            conditionalPanel("input['opt.downsampling.method']!='none'",
                                                                             sliderInput("downsampling", span("Down-sample Count",helpText("No display if more than four facets.")), 0, 100000, value=2000, step=1000)),
                                                            sliderInput("opt.expr.size", "Point Size of Maximum Expression", 1, 10, value=3, step=0.5),
                                                            checkboxInput("opt.scatter.gene.labels","Show Gene Labels on Scatter Plots", value=TRUE),
                                                            selectInput("opt.components", "Show Components", choices=c("Real"='real','Used for Clustering'='clustering','All'='all')),
                                                            fluidRow(column(7,checkboxInput("opt.use.cache", "Use Cached Plot Images", value=TRUE)), 
                                                                     column(2,actionButton("clear.cache","Clear Cache")))
                                                   )
                                       ),
                                       hr(),
                                       actionButton("dump","Save Debug State")
                          ),
                          
                          # Show a plot of the generated distribution
                          mainPanel(width=9, id="mainpanel",
                                    tabsetPanel(type="tabs",  
                                                tabPanel(span("Global Clusters"),
                                                         tabsetPanel(type="pills",
                                                                     tabPanel("tSNE",
                                                                              fluidRow(div(id="global-expression", class="scroll-area", 
                                                                                           plotDownload("tsne.global.cluster.label.dl"),
                                                                                           imageOutput("tsne.global.cluster.label", height="500px")))),
                                                                     tabPanel("Scatter", 
                                                                              fluidRow(div(id="global-scatter", class="scroll-area",
                                                                                           plotDownload("gene.expr.scatter.cluster.dl"),
                                                                                           imageOutput("gene.expr.scatter.cluster", height=500)))),
                                                                     tabPanel("Table",
                                                                              p("For clarity/debug. Probably not in final product."),
                                                                              fluidRow(div(id="dt-clusters", class="scroll-area", 
                                                                                           DT::dataTableOutput("dt.clusters"))))),
                                                         hr(),
                                                         tabsetPanel(type="pills",
                                                                     tabPanel("Differentially Over-Expressed Genes",
                                                                              fluidRow(class="table-area",
                                                                                       tableDownload("dt.cluster.markers.dl"),
                                                                                       column(2,p(style='height:1em'),actionButton("prev.cluster","<<<"), actionButton("next.cluster",">>>")),
                                                                                       column(4,uiOutput("current.cluster")),
                                                                                       column(4,offset=1, uiOutput("comparison.cluster"))),
                                                                              fluidRow(column(11,DT::dataTableOutput("dt.cluster.markers")))))),
                                                tabPanel("Subclusters",
                                                         tabsetPanel(type="pills",
                                                                     tabPanel("tSNE",
                                                                              fluidRow(div(id="local-expression", class="scroll-area",
                                                                                           conditionalPanel("input.showSubclustersInGlobal",
                                                                                                            plotDownload("tsne.global.subcluster.label.dl"),
                                                                                                            imageOutput("tsne.global.subcluster.label", height=500)),
                                                                                           conditionalPanel("!input.showSubclustersInGlobal",
                                                                                                            plotDownload("tsne.local.label.dl"),
                                                                                                            imageOutput("tsne.local.label", height=500)),
                                                                                           checkboxInput("showSubclustersInGlobal","Show Subclusters in Global Plot", value=FALSE)))),
                                                                     tabPanel("Scatter",
                                                                              fluidRow(div(id="global-scatter", class="scroll-area",
                                                                                           plotDownload("gene.expr.scatter.subcluster.dl"),
                                                                                           plotOutput("gene.expr.scatter.subcluster", height=500)))),
                                                                     tabPanel("Table",  
                                                                              fluidRow(div(id="dt-subclusters", class="scroll-area", DT::dataTableOutput("dt.subclusters"))))),
                                                         hr(),
                                                         tabsetPanel(type="pills",
                                                                     tabPanel("Differentially Over-Expressed Genes",
                                                                              fluidRow(class="table-area",
                                                                                       tableDownload("dt.subcluster.markers.dl"),
                                                                                       column(2,p(style='height:1em'),actionButton("prev.subcluster","<<<"), actionButton("next.subcluster",">>>")),
                                                                                       column(4,uiOutput("current.subcluster")),
                                                                                       column(4,offset=1, uiOutput("comparison.subcluster"))),
                                                                              fluidRow(column(11,DT::dataTableOutput("dt.subcluster.markers")))),
                                                                     tabPanel("Independent Components",
                                                                              fluidRow(DT::dataTableOutput("dt.components")))))
                                                # ,
                                                # tabPanel("Metacells")
                                    )
                          )
                          #    tabPanel("Community Annotations", p("A wiki-like editable annotation of cell type"))
                        )
               ),
               tabPanel("Data",
                        p("Links to data sets for downloading")),
               tabPanel("Help",
                        p("Explanations about using the app. We should also have some help info within the app panels.")),
               tabPanel("About",
                        p("Credits, manuscript, contact go here"))
    ),
    tags$script(src="http://code.jquery.com/ui/1.12.1/jquery-ui.min.js"),
    tags$script("jQuery(function (){ $('.scroll-area').resizable(); });")
  )
)
