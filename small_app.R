library(shiny)
library(tidyverse)
library(colourpicker) # you might need to install this.
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(SummarizedExperiment)
library(DT)
library(ggbeeswarm)

ui <- fluidPage(
  titlePanel("FINAL PROJECT"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput(inputId = "fgsea_file", label = "Load a CSV file", accept = ".csv")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(title = "Pathway table",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                sliderInput(inputId = "path_slid",label = "filter table by padj value", min = -48, max = 0, value = -20, step = 1),
                                radioButtons(inputId = "all_path", label = "select pathways",choices = c("All","Postive","Negative"), selected = NULL),
                                downloadButton(outputId = "download_fgsea_table", label = "Download")
                   ),
                   # Show a plot of the fgsea bars of top results
                   mainPanel(
                     div(DT::dataTableOutput("fgsea_filt_table"), style = "font-size:80%; width: 30%;")
                   )
                 )
        ),
        tabPanel(title = "Scatter Plot",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                sliderInput(inputId = "scatter_slid",label = "filter the plot by padj value", min = -48, max = 0, value = -20, step = 1)
                   ),
                   # Show a plot of the fgsea bars of top results
                   mainPanel( 
                     plotOutput("NES_scatter",width = "130%", height = "500px")
                   )
                 )
        )
      )
      
    )
  )
)


server <- function(input, output) {
  
  fgsea_data <- reactive({
    req(input$fgsea_file)
    data <-read.csv(input$fgsea_file$datapath, header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()
    return(data)
  })
  
  #Generate the filter pathway tables in fgsea
  gsea_table <- reactive({
    filter_df<-fgsea_data()
    filtered_df <- dplyr::filter(filter_df, padj < 10^(input$path_slid))
    if (input$all_path == "All") {
      filtered_df <- filtered_df
    } else if (input$all_path == "Positive") {
      filtered_df <- dplyr::filter(filtered_df, NES > 0)
    } else if (input$all_path == "Negative") {
      filtered_df <- dplyr::filter(filtered_df, NES < 0)
    }
    filtered_df <- filtered_df %>%
      dplyr::mutate(pval = formatC(ifelse(is.na(pval), 0, pval), format = "e"),
                    padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
    return(filtered_df)
  })
  
  #Output download fgsea table
  output$download_fgsea_table <- downloadHandler(
    filename = function() {
      paste("fgsea_table", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(gsea_table(), file, row.names = FALSE)
    }
  )
  
  output$fgsea_filt_table <- DT::renderDataTable({
    filter_df <- fgsea_data()
    filtered_df <- dplyr::filter(filter_df, padj < 10^(input$path_slid))
    if (input$all_path == "All") {
      filtered_df <- filtered_df
    } else if (input$all_path == "Positive") {
      filtered_df <- dplyr::filter(filtered_df, NES > 0)
    } else if (input$all_path == "Negative") {
      filtered_df <- dplyr::filter(filtered_df, NES < 0)
    }
    filtered_df <- filtered_df %>%
      dplyr::mutate(pval = formatC(ifelse(is.na(pval), 0, pval), format = "e"),
                    padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
    return(DT::datatable(filtered_df))
  })
    
  output$NES_scatter <- renderPlot({
    dataf <- fgsea_data()
    slider<-input$scatter_slid
    df <- dplyr::mutate(dataf, new_padj = -log10(padj)) %>%
      dplyr::mutate(status = dplyr::case_when(padj < 10^(slider) ~ "TRUE",
                                              padj >= 10^(slider) ~ "FALSE"))
    # specify color based on the slider value
    df$colors <- ifelse(df$status == "FALSE", "orange", "grey")
    # plotting scatter plot
    scatter <- ggplot(df, aes(x = NES, y = new_padj, color = colors)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("orange","grey"),
                         labels = c("TRUE", "FALSE")) +
      labs(x = "NES", y = "-log10(padj)",color = paste0( "padj < 10^",slider )) +
      theme_bw()+
      theme(legend.position = "bottom") # move legend to bottom of plot
    return(scatter)
  }) 
  
  
}

shinyApp(ui, server)
