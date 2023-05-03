# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(tidyverse)
library(colourpicker) # you might need to install this.
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(SummarizedExperiment)
library(DT)

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
    ),
    tabPanel(title = "DE",
       sidebarLayout(
         sidebarPanel(
           #input count matrix
           fileInput(inputId = "deseq_file", label = "Load a CSV file", accept = ".csv"),
           # Add submit button
           submitButton(text = "Submit",icon = icon("chart-line"))
         ),
         # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "DE input file",
              DT::dataTableOutput("DE_summary")
            ),
            tabPanel(title = "DE result",
              sidebarPanel(
               radioButtons(inputId = "x_axis", label = "select x variable", choices = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "log2FoldChange"), #NULL ),
               radioButtons(inputId = "y_axis", label = "select y variable", choices = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "padj"),#NULL),
               # Add color inputs
               colourInput(inputId = "base", label = "Choose color 1", value = "#07B377"),
               colourInput(inputId = "highlight", label = "Choose color 2", value = "#F5E149"),
               # Add slider inputs
               sliderInput(inputId = "padj_slider",label = "Choose a padj value as threshold", min = -35, max = 0, value = -8, step = 1),
               #Add a submit buttom
               submitButton(text = "plot",icon = icon("chart-line")) #style = ......
              ),
              mainPanel(
               tabsetPanel(
                 tabPanel(title = "Volcano plot",
                          plotOutput("volcano",width = "90%", height = "550px")
                 ),
                 tabPanel(title = "Padj filtered table",
                          DT::dataTableOutput("volcano_table")
                 )
               )
              )
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
          fileInput(inputId = "fgsea_file", label = "Load a CSV file", accept = ".csv"),
          # Add submit button
          submitButton(text = "Submit",icon = icon("chart-line"))
        ),
       # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "Pathway barplot",
              sidebarLayout(
                sidebarPanel( 
                  sliderInput(inputId = "pth_threshold",label = "Top results by padj value", min = -48, max = 0, value = -20, step = 1),
                  submitButton(text = "Submit",icon = icon("chart-line"))
                ),
                # Show a plot of the fgsea bars of top results
                mainPanel( 
                  plotOutput("fgsea_bars")
                )
              )
            ),
            tabPanel(title = "Pathway table",
              sidebarLayout(
                sidebarPanel(sliderInput(inputId = "path_slid",label = "filter table by padj value", min = -48, max = 0, value = -20, step = 1),
                            radioButtons(inputId = "all_path", label = "select pathways",choices = c("All","Postive","Negative"), selected = NULL), 
                            downloadButton(outputId = "download_fgsea_table", label = "Download")
                ),
                # Show a plot of the fgsea bars of top results
                mainPanel(
                  DT::dataTableOutput("fgsea_filt_table")
                )
              )
            ),
            tabPanel(title = "Scatter Plot",
              sidebarLayout(
                sidebarPanel(sliderInput(inputId = "scatter_slid",label = "filter the plot by padj value", min = -48, max = 0, value = -20, step = 1)
                ),
                # Show a plot of the fgsea bars of top results
                mainPanel( 
                  plotOutput("NES_scatter")
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
  options(shiny.maxRequestSize = 30*1024^2) # set max file size to 30 MB
  #Read in sample data
  sample_data <- reactive({
    req(input$sample_file)
    data <-read.csv(input$sample_file$datapath, sep="\t", header = FALSE, stringsAsFactors = FALSE)%>%as_tibble()
    data_info <- data %>% 
      mutate(V1 = apply(data["V1"],1, function(x) gsub("!", "", x))) %>% 
      as.data.frame()
    data_info<-data_info[c(1:3,6,8:14),]%>%
      mutate(Columns = c("GEO accession","Status","Submission date", "Channel count","Organism_source","Tissue_source","Diagnosis", "Post_mortem_interval","Age_of_death","RNA_integrity_number","mRNA_seq_reads"), .before = 1) %>%
      select(-V1)%>%
      t()
    colnames(data_info) <- data_info[1,]
    data_info<- data_info[-1,]%>%
      apply(2, function(x) gsub("tissue: ", "", x)) %>%
      apply(2, function(x) gsub("diagnosis: ", "", x)) %>%
      apply(2, function(x) gsub("pmi: ", "", x)) %>%
      apply(2, function(x) gsub("age of death: ", "", x)) %>%
      apply(2, function(x) gsub("rin: ", "", x)) %>%
      apply(2, function(x) gsub("mrna-seq reads: ", "", x)) %>%
      as.data.frame() 
    #names(unknowndata_info) <- as.character(unknowndata_info[1,])
    data_info <-mutate(data_info,across(c(4,8:11), as.double))
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
  #generate PCA beeswarmplot
  plot_beeswarm <-function(data){
    pca_results <- prcomp(scale(t(data[,-1]%>%as.data.frame())), center=FALSE, scale=FALSE)
    plot_tibble <- as_tibble(pca_results$x)%>%
      add_column( sample = rownames(pca_results$x), .after=0)
    meta <- tibble(sample = rownames(pca_results$x)) %>%
      mutate(Diagnosis = if_else(row_number() <= 49, "normal", "Huntington's Disease"))
    biplot <- dplyr::left_join(plot_tibble, meta, by= "sample")
    biplot$Diagnosis <- factor(ifelse(biplot$Diagnosis == "normal" & seq_len(nrow(biplot)) <= 49, "normal", "Huntington's Disease"))
    biplot_select <- dplyr::select(biplot, 1:21, 71)
    # Define a custom color palette with repeated colors
    my_colors <- rep(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5"), 2)
    beeswarm_plot <- biplot_select %>%
      gather(key = "PC", value = "value", PC1:PC20) %>%
      ggplot(aes(x = factor(PC, levels = paste0("PC", 1:20)), y = value, color = PC)) +
      geom_quasirandom(size = 0.6, width = .3) +
      scale_color_manual(values = my_colors) +
      theme_classic()+
      labs(x = "PCs", y = "Values")
    beeswarm_plot
  }
  #Read in DESeq data
  deseq_data <- reactive({
    req(input$deseq_file)
    data <-read.table(input$deseq_file$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()%>%dplyr::rename(gene= X)
    return(data)
  })
  
  #Generate Volcano plot 
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    # modify the dataframe
    #y_name <- gsub('"','',y_name)
    df <- dplyr::mutate(dataf, new_y_name = -log10(!!sym(y_name))) %>%
      dplyr::mutate(status = dplyr::case_when(padj < 10^(slider) ~ "TRUE",
                                              padj >= 10^(slider) ~ "FALSE",
                                              TRUE ~ "NA"))
    
    # specify color based on the slider value
    df$colors <- ifelse(df$status == "FALSE", color1, color2)
    #plotting volcano plot
    volcano <- ggplot(df, aes(x = !!sym(x_name), y = new_y_name,color = colors)) + 
      geom_point(size = 1) +
      scale_color_manual(values = c(color1, color2,"grey"),
                         labels = c("FALSE", "TRUE", "NA")) +
      labs(x = x_name, y = paste0("-log10(", y_name, ")"),color = paste0( y_name, "< 10^",slider )) +
      theme_bw()+
      theme(legend.position = "bottom") # move legend to bottom of plot
    
    return(volcano)
  }
  
  #Generate padj filtered table in DE tab
  draw_table <- function(dataf, slider) {
    filtered_df <-dplyr::filter(dataf, padj < 10^(slider))
    filtered_df <- dplyr::rename(filtered_df)
    formatted_df <- dplyr::mutate(filtered_df, pvalue = formatC(ifelse(is.na(pvalue), 0, pvalue), format = "e"),
                                  padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
    return(formatted_df)
  }
  #Read in fgsea data
  fgsea_data <- reactive({
    req(input$fgsea_file)
    data <-read.csv(input$fgsea_file$datapath, header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()
    return(data)
  })
  
  #Generate barplot for top pathways of fgsea
  fgsea_top_pathways <- function(fgsea_results, threshold){
    top_positive_nes <- fgsea_results %>%
      dplyr::filter(padj < 10^(threshold) & NES > 0)
    top_negative_nes <- fgsea_results %>%
      dplyr::filter(padj < 10^(threshold) & NES <0)
    NES_barplot <- dplyr::bind_rows(top_positive_nes, top_negative_nes)%>%
      ggplot() +
      geom_col(aes(x=reorder(pathway,+NES), y=NES, fill = NES > 0))+
      scale_fill_manual(values =c('TRUE' = 'red', 'FALSE' = 'blue')) +
      theme_minimal() +
      coord_flip()+
      theme(legend.position = "none", axis.text.y = element_text( hjust =1 ,size= 3),axis.title.x = element_text(size = 6))+ 
      labs(title="fgsea results for Hallmark MSigDB gene sets", x= "",y= "Normalized Enrichment Score (NES)")
    return(NES_barplot)
  }
  
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
  
  #Generate scatter plot in GSEA
  scatter_plot <- function(dataf, slider) {
    # modify the dataframe
    df <- dplyr::mutate(dataf, padj = -log10(padj)) %>%
      dplyr::mutate(status = dplyr::case_when(padj < 10^(slider) ~ "TRUE",
                                              padj >= 10^(slider) ~ "FALSE",
                                              TRUE ~ "NA"))
    
    # specify color based on the slider value
    df$colors <- ifelse(df$status == "FALSE", "orange", "grey")
    # plotting volcano plot
    scatter <- ggplot(df, aes(x = NES, y = padj, color = colors)) +
      geom_point(size = 1) +
      scale_color_manual(values = c("orange","grey"),
                         labels = c("FALSE", "TRUE")) +
      labs(x = "NES", y = "-log10(padj)",color = paste0( "padj < 10^",slider )) +
      theme_bw()+
      theme(legend.position = "bottom") # move legend to bottom of plot
    
    return(scatter)
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
  
  #Output filtering count heatmap
  output$clus_heatmap <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    num_matrix <- filter_res(dataf, input$slid_var, input$slid_zero)
    plot_heatmap(num_matrix)
  }) 
  
  #Output count PCA beeswarmplot
  output$pca_plot <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    beeswarm_plot <- plot_beeswarm(dataf)
    beeswarm_plot
  }) 
  
  #output DESeq summary table
  output$DE_summary <- DT::renderDataTable({
    dataf <- deseq_data()
    if(is.null(dataf))
      return(null)
    dataf
  }) 
  
  #output DE volcano plot
  output$volcano <- renderPlot({
    dataf <- deseq_data()
    if(is.null(dataf))
      return(null)
    plot <- volcano_plot(dataf, input$x_axis, input$y_axis, input$padj_slider, input$base, input$highlight)
    #ggplotly(plot) # turn the static plot into interactive but correspoding to use use plotlyOutput() and renderPlotly()
    plot
  }) 
  
  # output the padj filtered table in DE tab
  output$volcano_table <- DT::renderDataTable({
    dataf <- deseq_data()
    if(is.null(dataf))
      return(null)
    result_tab <- draw_table(dataf, input$padj_slider)
    result_tab
  }) 
  
  #Output fgsea top barplot
  output$fgsea_bars <- renderPlot({
    dataf <- fgsea_data()
    if(is.null(dataf))
      return(null)
    bar_plot <- fgsea_top_pathways(dataf,input$pth_threshold)
    bar_plot
  }) 
  
  #output fgsea pathway filter table
  #let the input be reactive so when we move the bar, the output would change
  output$fgsea_filt_table <- DT::renderDataTable({
    gsea_table()
  }) 
  #Output fgsea filter scatter plot
  output$NES_scatter <- renderPlot({
    dataf <- fgsea_data()
    if(is.null(dataf))
      return(null)
    scatter_plot <- scatter_plot(dataf,input$scatter_slid)
    scatter_plot
  }) 
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
