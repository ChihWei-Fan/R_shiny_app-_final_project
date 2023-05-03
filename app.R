#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggplot2)
library(colourpicker)

# Define UI for application that draws a histogram
ui <- fluidPage(
  #Application title
  titlePanel("Old Faithful Geyser Data"),
  # 4 tabs
  tabsetPanel(
    tabPanel(title = "Sample",

     sidebarLayout(
       sidebarPanel(
         fileInput(inputId = "sample_file", label = "Load a CSV file", accept = ".csv"),
         submitButton(text = "Submit",icon = icon("chart-line"))
       ),
       
       # Show a plot of the generated distribution
       mainPanel(
         tabsetPanel(
           tabPanel(title = "Summary",
                    tableOutput("summary_file")
           ),
           tabPanel(title = "Table",
                    tableOutput("summary_table")
           ),
           tabPanel(title = "Plot",
                    plotOutput("plotname",width = "90%", height = "550px")
           )
         )
       )
     )
    ),
    tabPanel(title = "Counts",
      sidebarLayout(
       sidebarPanel(
         #input count matrix
         fileInput(inputId = "count_file", label = "Load a CSV file", accept = ".csv"),
         # Add slider inputs
         sliderInput(inputId = "slid_var",label = "Choose a threshold value to include gene at least X percentile of variance", min = -300, max = 0, value = -150, step = 1),
         sliderInput(inputId = "slid_zero",label = "Choose a threshold value to include gene in X sample are non-zero", min = -300, max = 0, value = -150, step = 1),
         submitButton(text = "Submit",icon = icon("chart-line"))
       ),
       
       # Show a plot of the generated distribution
       mainPanel(
         tabsetPanel(
           tabPanel(title = "Filter effect",
                    tableOutput("summary_file")
           ),
           tabPanel(title = "Scatter plot",
                    plotOutput("count_scatter")
           ),
           tabPanel(title = "Heatmap",
                    plotOutput("clus_heatmap",width = "90%", height = "550px")
           ),
           tabPanel(title = "PCA",
                    plotOutput("pca_plot",width = "90%", height = "550px")
           )
         )
       )
      )
    ),
    tabPanel(title = "DE",
     sidebarLayout(
       sidebarPanel(
         #input count matrix
         fileInput(inputId = "file", label = "Load a CSV file", accept = ".csv"),
         # Add submit button
         submitButton(text = "Submit",icon = icon("chart-line"))
       ),
       
       # Show a plot of the generated distribution
       mainPanel(
         tabsetPanel(
           tabPanel(title = "DE result",
                    tableOutput("summary_file")
           ),
           tabPanel(title = "NOT KNOWN",
                    plotOutput("diag_scatter")
           )
         )
       )
     )
    ),
    tabPanel(title = "GSEA",
      # Use DGE results to compute gene set enrichment analysis with fgsea
      sidebarLayout(
        sidebarPanel(
          #input count matrix
          fileInput(inputId = "file", label = "Load a CSV file", accept = ".csv"),
          # Add submit button
          submitButton(text = "Submit",icon = icon("chart-line"))
        ),
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "TOP result",
             sidebarLayout(
               sidebarPanel( sliderInput(inputId = "slid_var",label = "Top results by padj value", min = 0, max = 100, value = 60, step = 1)
               ),
               # Show a plot of the fgsea bars of top results
               mainPanel( plotOutput("fgsea_bars")
               )
             )
            ),
            tabPanel(title = "Table",
              sidebarLayout(
                sidebarPanel(sliderInput(inputId = "slid_var",label = "filter table by padj value", min = 0, max = 100, value = 60, step = 1),
                             radioButtons(inputId = "all_pth", label = "select pathways",choices = c("All","Postive","Negative"), selected = NULL), 
                             downloadButton(outputId = "download_fgses_table", label = "Download")
                ),
               # Show a plot of the fgsea bars of top results
                mainPanel( tableOutput("fgsea_sort_table")
                )
              )
            ),
            tabPanel(title = "Plots",
              sidebarLayout(
               sidebarPanel(sliderInput(inputId = "slid_var",label = "filter table by padj value", min = 0, max = 100, value = 60, step = 1)
               ),
               # Show a plot of the fgsea bars of top results
               mainPanel( plotOutput("NES_scatter")
               )
              )
            )
          )
        )
      )
    )
  )  
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  plot_variance_vs_median <- function(data, scale_y_axis=FALSE) {
    new_tib <-tibble(count_median = apply(data,1,median))%>%
      add_column(variance = apply(data, 1, var), rank = rank(.))%>%
      ggplot(aes(x= rank, y= variance))+
      geom_point()+
      geom_smooth()+
      labs(title=title,x="Rank(Median)", y = "Variance")+
      scale_y_log10()
    
    return(new_tib)
  }
  
  sample_data <- reactive({
    req(input$sample_file)
    data <- read_data(input$sample_file$datapath)
    return(data)
  })
  
  count_data <- reactive({
    req(input$count_file)
    data <-read.table(input$count_file$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()%>%
      rename(gene = X)
    return(data)
  })
  
  output$count_scatter <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    plot <- plot_variance_vs_median(dataf)
    #ggplotly(plot) # turn the static plot into interactive but correspoding to use use plotlyOutput() and renderPlotly()
    plot
  }) 
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
