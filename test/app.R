# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(tidyverse)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(gridExtra)
library(RColorBrewer)
library(gplots)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 final Project"),
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
                            tableOutput("sample_summary")
                   ),
                   tabPanel(title = "Table",
                            tableOutput("sample_table")
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
           sliderInput(inputId = "slid_var",label = "Choose a threshold value to include gene at least X percentile of variance", min = 0, max = 100, value = 10, step = 1),
           sliderInput(inputId = "slid_zero",label = "Choose a threshold value to include gene in X sample are non-zero", min = 0, max = 24, value = 5, step = 1),
           submitButton(text = "Submit",icon = icon("chart-line"))
         ),
         
         # Show a plot of the generated distribution
         mainPanel(
           tabsetPanel(
             tabPanel(title = "Filter effect",
                      tableOutput("filter_count")
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
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 30*1024^2) # set max file size to 30 MB
  #Read in sample data
  sample_data <- reactive({
    req(input$sample_file)
    data <- read.csv(input$sample_file$datapath, sep="\t", header=FALSE, stringsAsFactors=FALSE) %>%
      as_tibble()
    data <- data %>% mutate(V1 = apply(data["V1"],1, function(x) gsub("!", "", x)))%>%
      slice(1:3,6,8:14) %>%
      mutate(Columns = c("GEO accession","Status","Submission date", "Channel count","Organism_source","Tissue_source","Diagnosis", "Post_mortem_interval","Age_of_death","RNA_integrity_number","mRNA_seq_reads"), .before = 1) %>%
      select(-V1) %>%
      t() %>%
      apply(2, function(x) gsub("tissue: ", "", x)) %>%
      apply(2, function(x) gsub("diagnosis: ", "", x)) %>%
      apply(2, function(x) gsub("pmi: ", "", x)) %>%
      apply(2, function(x) gsub("age of death: ", "", x)) %>%
      apply(2, function(x) gsub("rin: ", "", x)) %>%
      apply(2, function(x) gsub("mrna-seq reads: ", "", x)) %>%
      as.data.frame() 
    names(data) <- as.character(data[1,])
    data <- data[-1,] %>%
      mutate(across(c(4,8:11), as.double))
  })
  
  summary_tablef <- function(data) {
    # Count number of rows and columns
    n_rows <- nrow(data)
    n_cols <- ncol(data)
    
    # Create a data frame with column information
    col_info <- data.frame(
      Column_name = names(data),
      Data_type = sapply(data, class),
      stringsAsFactors = FALSE
    )
    # Add mean or distinct values for columns
    for (i in 1:n_cols) {
      if (is.numeric(data[[i]])) {
        col_info$mean[i] <- mean(data[[i]], na.rm = TRUE)
      } 
      else { col_info$mean[i] <- "N/A"
      }
    }
    return(col_info%>%as_tibble())
  }
  
  #generate sample summary table
  output$sample_summary <- renderTable({
    dataf <- sample_data()
    if(is.null(dataf))
      return(null)
    samplesum_tab <- summary_tablef(dataf)
    samplesum_tab
  }) 
  
  #Generate sample file as table
  output$sample_table <- renderTable({
    dataf <- sample_data()
    if(is.null(dataf))
      return(null)
    dataf
  }) 
  
  
  #Read in count data
  count_data <- reactive({
    req(input$count_file)
    data <-read.table(input$count_file$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()%>%
      rename(gene = X)
    return(data)
  })
  
  filter_table <- function(data, pass_filter2) {
    count_sum <-tibble(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0))) %>%
      mutate(determine = ifelse(num_zeros >= pass_filter2, "pass", "fail")) %>%
      summarise(num_passing = sum(determine == "pass"),
                perc_passing = round(num_passing/nrow(.) * 100, 2))
    # Create a data frame with filter information
    col_info <- data.frame(
      "Number of samples" = ncol(data),
      "Number of genes" = nrow(data),
      "Number and % of genes passing current filter" = paste0("Num of genes ",count_sum[[1]], ", percentage of genes",count_sum[[2]],"%"),
      "Number and % of genes not passing current filter"= paste0("Num of genes ",nrow(data)-count_sum[[1]], ", percentage of genes",100-count_sum[[2]],"%"),
      stringsAsFactors = FALSE
    )
    return(col_info)
  }
  
  plot_variance_vs_median <- function(data, pass_filter, scale_y_axis=FALSE) {
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(variance = apply(data[,-1], 1, var), rank = rank(.)) %>%
      mutate(determine = ifelse(variance >= pass_filter, "pass", "fail")) %>%
      ggplot(aes(x= rank, y= variance, color=determine)) +
      geom_point() +
      scale_color_manual(values = c("fail" = "lightblue", "pass" = "darkblue")) +
      geom_smooth() +
      labs(title= "Median Count vs Variance", x="Rank(Median)", y = "Variance", color="Filter") +
      scale_y_log10()
    return(new_tib)
  }
  
  plot_variance_vs_nonzero <- function(data, pass_filter2, scale_y_axis=FALSE) {
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0)), rank = rank(.)) %>%
      mutate(determine = ifelse(num_zeros >= pass_filter2, "pass", "fail")) %>%
      ggplot(aes(x= rank, y= num_zeros, color=determine)) +
      geom_point() +
      scale_color_manual(values = c("fail" = "lightblue", "pass" = "darkblue")) +
      geom_smooth() +
      labs(title= "Median Count vs Number of zeros", x="Rank(Median)", y = "Number of zeros", color="Filter")
    return(new_tib)
  }
  #Fitler the unqualified genes
  filter_res <- function(data, pass_filter1, pass_filter2) {
    filt_res <-tibble(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0)),
                      variance = apply(data[,-1], 1, var)) %>%
      mutate(determine_var = ifelse(variance >= pass_filter1, "pass", "fail"),
             determine_0 = ifelse(num_zeros >= pass_filter2, "pass", "fail"))%>%
      filter(determine_var == "pass" & determine_0 == "pass")
    return (filt_res)
  }
  #generate filter matrix for heatmap
  filter_res <- function(data, pass_filter1, pass_filter2) {
    filt_res <- data %>% 
      mutate(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0)),
             variance = apply(data[,-1], 1, var)) %>%
      filter(variance > pass_filter1 & num_zeros > pass_filter1)%>% as.data.frame()
    rownames(filt_res) <-filt_res[,1]
    return (filt_res)
  }
  #generate count heatmap after filtering
  plot_heatmap <- function(filter_data) {
    coul <- rev(brewer.pal(11, 'RdBu'))
    num_matrix <- filter_data[-c(1, 71, 72)] %>% as.matrix()%>% log2()
    num_matrix[!is.finite(num_matrix)] <- NA
    heatmap.2(num_matrix, col = coul,trace = "none")
  }
  
  #output count filter table
  output$filter_count <- renderTable({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    result_tab <- filter_table(dataf, input$slid_zero)
    result_tab
  }) 
  
  #output count scatter plot
  output$count_scatter <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    plot1 <- plot_variance_vs_median(dataf,input$slid_var)
    plot2 <- plot_variance_vs_nonzero (dataf,input$slid_zero)
    #ggplotly(plot) # turn the static plot into interactive but correspoding to use use plotlyOutput() and renderPlotly()
    grid.arrange(plot1, plot2, nrow=2)
  }) 
  
  #
  output$clus_heatmap <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    num_matrix <- filter_res(dataf, input$slid_var, input$slid_zero)
    plot_heatmap(num_matrix)
  }) 
  
  
  
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
