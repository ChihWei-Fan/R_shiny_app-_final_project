# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(tidyverse)
library(colourpicker) # you might need to install this.
library(gridExtra)
library(RColorBrewer)
library(gplots)
library(SummarizedExperiment)
library(DT)
library(ggbeeswarm)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Final Project"),
  h4("Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"),
  tabsetPanel(
    tabPanel(title = "Sample",
      sidebarLayout(
       sidebarPanel( width = 3,
         fileInput(inputId = "sample_file", label = "Load a CSV file", accept = ".csv")
       ),
               
               # Show a plot of the generated distribution
        mainPanel(
          
          tabsetPanel(
            tabPanel(title = "Summary of Sample File",
                    tableOutput("sample_summary")
            ),
            tabPanel(title = "Sample Content Table",
              div(DT::dataTableOutput("sample_table"), style = "font-size:80%; width: 30%;")
            ),
            tabPanel(title = "Sample Distribution Plot",
              sidebarPanel( width = 3,
                radioButtons(inputId = "sample_x", label = "select y variable", choices = c("Post mortem interval","Age of death","RNA integrity number","mRNA-seq reads"), selected = "Age of death"),
                  submitButton(text = "Submit",icon = icon("print"))
              ),
              mainPanel(
                plotOutput("sample_plot")
              )
            )
          )
        )
      )
    ),  
    
    tabPanel(title = "Counts Matrix Exploration",
      sidebarLayout(
        sidebarPanel( width = 3,
           #input count matrix
           fileInput(inputId = "count_file", label = "Load a CSV file", accept = ".csv"),
           # Add slider inputs
           sliderInput(inputId = "slid_var",label = "Choose a threshold value to include genes that have at least X percentile of variance", min = 1, max = 100, value = 80, step = 1),
           sliderInput(inputId = "slid_zero",label = "Choose a threshold value to include genes that have at least X samples with non-zero values", min = 0, max = 69, value = 60, step = 1),
           submitButton(text = "Submit",icon = icon("file"))
         ),
         
         # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "Filter Results",
              tableOutput("filter_count")
            ),
            tabPanel(title = "Scatter Plot of Fitered Count Matrix",
              plotOutput("count_scatter")
            ),
            tabPanel(title = "Heatmap of Fitered Count Matrix",
              plotOutput("clus_heatmap",width = "80%", height = "500px")
            ),
            tabPanel(title = "PCA of Count Matrix",
              sidebarPanel(width = 3,
                sliderInput("top_PC",label = "Select the TOP PCs you want to plot", min = 3, max = 45, value =8, step = 1),
                submitButton(text = "Submit",icon = icon("chart-line"))
              ),
              mainPanel(
                plotOutput("pca_plot",width = "120%", height = "400px")
              )
            )
          )
        )
      )
    ),
    tabPanel(title = "Differential Expression Analysis",
       sidebarLayout(
         sidebarPanel(width = 3,
           #input count matrix
           fileInput(inputId = "deseq_file", label = "Load a CSV file", accept = ".csv"),
         ),
         # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "DE Input File",
              div(DT::dataTableOutput("DE_summary"), style = "font-size:80%; width: 30%;")
            ),
            tabPanel(title = "DE Results",
              sidebarPanel(width = 3,
               radioButtons(inputId = "x_axis", label = "Select a variable on X axis", choices = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "log2FoldChange"), #NULL ),
               radioButtons(inputId = "y_axis", label = "Select a variable on Y axis", choices = c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"), selected = "padj"),#NULL),
               # Add color inputs
               colourInput(inputId = "base", label = "Choose color 1", value = "#07B377"),
               colourInput(inputId = "highlight", label = "Choose color 2", value = "#F5E149"),
               # Add slider inputs
               sliderInput(inputId = "padj_slider",label = "Choose a padj value as threshold", min = -35, max = 0, value = -8, step = 1),
               #Add a submit buttom
               submitButton(text = "plot",icon = icon("folder")) 
              ),
              mainPanel(
                tabsetPanel(
                  tabPanel(title = "Volcano plot",
                    plotOutput("volcano",width = "110%", height = "500px")
                  ),
                  tabPanel(title = "Padj filtered table",
                    div(DT::dataTableOutput("volcano_table"), style = "font-size:80%; width: 30%;")
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
        sidebarPanel(width = 3,
          #input count matrix
          fileInput(inputId = "fgsea_file", label = "Load a CSV file", accept = ".csv")
        ),
       # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel(title = "Pathway Barplot",
              sidebarLayout(
                sidebarPanel( width = 3,
                  sliderInput(inputId = "pth_threshold",label = "Top results by padj value smaller than 10^(X)", min = -48, max = 0, value = -20, step = 1),
                  submitButton(text = "Submit",icon = icon("bars"))
                ),
                # Show a plot of the fgsea bars of top results
                mainPanel( 
                  plotOutput("fgsea_bars",width = "130%", height = "500px")
                )
              )
            ),
            tabPanel(title = "Pathway Table",
                     sidebarLayout(
                       sidebarPanel(width = 3,
                                    sliderInput(inputId = "path_slid",label = "Filtered table by padj value", min = -48, max = 0, value = -20, step = 1),
                                    radioButtons(inputId = "all_path", label = "Select pathways",choices = c("All","Postive","Negative"), selected = NULL),
                                    submitButton(text = "Submit",icon = icon("filter")),
                                    downloadButton(outputId = "download_fgsea_table", label = "Download")
                       ),
                       # Show a plot of the fgsea bars of top results
                       mainPanel(
                         div(DT::dataTableOutput("fgsea_filt_table"), style = "font-size:80%; width: 30%;")
                       )
                     )
            ),
            tabPanel(title = "Pathway Scatter Plot",
                     sidebarLayout(
                       sidebarPanel(width = 3,
                                    sliderInput(inputId = "scatter_slid",label = "filter the plot by padj value", min = -48, max = 0, value = -20, step = 1),
                                    submitButton(text = "Submit",icon = icon("refresh"))
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
      mutate(Columns = c("GEO accession","Status","Submission date", "Channel count","Organism source","Tissue source","Diagnosis", "Post mortem interval","Age of death","RNA integrity number","mRNA-seq reads"), .before = 1) %>%
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
    rownames(data_info) <- NULL
    data_info <-mutate(data_info,across(c(4,8:11), as.double))
  })
  
  summary_tablef <- function(data) {
    # Count number of rows and columns
    n_rows <- nrow(data[,-1])
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
  
  #Generate histogram of sample
  histogram <- function(df, col_name) {
    binwidth <- if_else(!!sym(col_name) == "mRNA-seq reads", 20000000, 2)
    ggplot(df, aes(x = !!sym(col_name))) +
      geom_histogram(binwidth = binwidth, color = "black", fill = "purple") +
      theme_linedraw()+
      labs(x = col_name, y = "Count")
  }
  
  #Read in count data
  count_data <- reactive({
    req(input$count_file)
    data <-read.table(input$count_file$datapath, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
      as_tibble()%>%dplyr::rename(gene = X)
    return(data)
  })
  
  #Count data filtering table
  filter_table <- function(data, pass_filter1, pass_filter2) {
    count_sum <-data[,-1]
    count_sum$non_zeros <- apply(data[,-1], 1, function(x) sum(x != 0))
    geneVar <- apply(data[,-1], 1, var)
    #Count the percentile variance
    xPercentile <- quantile(geneVar, pass_filter1/100) #assume the filter is 1 to 100
    #find out which rows pass the filter1
    pass_xPercentile <- which(geneVar >= xPercentile)
    #subset the matrix to let the rows that pass remain
    count_sum <- count_sum[pass_xPercentile, ]
    count_sum$pass_zero <- ifelse(count_sum$non_zeros >= pass_filter2, "Pass", "Fail")
    num_passing <- sum(count_sum$pass_zero == "Pass")
    perc_passing <- round(num_passing/nrow(data[,-1])* 100, 2)
    # Create a data frame with filter information
    col_info <- tibble(
      metric = c("Number of samples", "Number of genes ", "Number of genes passing", " % of genes passing", "Number of genes not passing","% of genes not passing"),
      value = c(ncol(data[,-1]), nrow(data[,-1]) , num_passing, paste0(perc_passing,"%"), nrow(data[,-1])-num_passing, paste0(100-perc_passing,"%"))
    )
    return(col_info)
  }
  
  plot_variance_vs_median <- function(data, pass_filter1, scale_y_axis=FALSE) {
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(variance = apply(data[,-1], 1, var), rank = rank(.))
    #Count the percentile variance
    xPercentile <- quantile(new_tib$variance, pass_filter1/100)
    new_tib <- new_tib%>%  
      mutate(determine = ifelse(variance >= xPercentile, "Pass", "Fail")) %>%
      ggplot(aes(x= rank, y= variance, color=determine)) +
      geom_point() +
      scale_color_manual(values = c("Fail" = "lightblue", "Pass" = "darkblue")) +
      geom_smooth() +
      theme_classic()+
      labs(title= "Median Count vs Variance", x="Rank(Median)", y = "Variance", color="Filter") +
      scale_y_log10()
    return(new_tib)
  }
  
  plot_nonzero_vs_median <- function(data, pass_filter2, scale_y_axis=FALSE) {
    new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
      add_column(non_zeros = apply(data[, -1], 1, function(x) sum(x != 0)), rank = rank(.)) %>%
      mutate(determine = ifelse(non_zeros >= pass_filter2, "Pass", "Fail"),
             num_zeros = ncol(data)-non_zeros) %>%
      ggplot(aes(x= rank, y= num_zeros, color=determine)) +
      geom_point() +
      scale_color_manual(values = c("Fail" = "lightblue", "Pass" = "darkblue")) +
      geom_smooth() +
      theme_classic()+
      labs(title= "Median Count vs Number of zeros", x="Rank(Median)", y = "Number of zeros", color="Filter")
    return(new_tib)
  }

  #generate filter matrix for heatmap
  filter_res <- function(data, pass_filter1, pass_filter2) {
    filt_res <- data %>% 
      mutate(non_zeros = apply(data[, -1], 1, function(x) sum(x != 0)),
             variance = apply(data[,-1], 1, var))
    xPercentile <- quantile(filt_res$variance, pass_filter1/100)
    filt_res <- filt_res %>% filter(variance > xPercentile & non_zeros > pass_filter2)%>% as.data.frame()
    rownames(filt_res) <- filt_res[,1]
    return (filt_res[,-c(1,71,72)])
  }
  
  #generate count heatmap after filtering
  plot_heatmap <- function(filter_data) {
    coul <- rev(brewer.pal(11, 'RdBu'))
    num_matrix <- filter_data %>% as.matrix() %>% log2()
    num_matrix[!is.finite(num_matrix)] <- NA
    heatmap.2(num_matrix, col = coul, trace = "none",xlab = "Samples", ylab = "Genes",margins = c(5, 8),key = TRUE, key.title = "Expression level", key.xlab = "Expression", key.ylab = NULL)
  }
  #generate PCA beeswarmplot
  plot_beeswarm <- function(data, N) {
    pca_results <- prcomp(scale(t(data[,-1]%>%as.data.frame())), center=FALSE, scale=FALSE)
    plot_tibble <- as_tibble(pca_results$x) %>%
      add_column(sample = rownames(pca_results$x), .after = 0)
    meta <- tibble(sample = rownames(pca_results$x)) %>%
      mutate(Diagnosis = if_else(row_number() <= 49, "normal", "Huntington's Disease"))
    biplot <- dplyr::left_join(plot_tibble, meta, by = "sample")
    biplot$Diagnosis <- factor(ifelse(biplot$Diagnosis == "normal" & seq_len(nrow(biplot)) <= 49, "normal", "Huntington's Disease"))
    biplot_select <- dplyr::select(biplot, 1:N+1, 71)
    # Define a custom color palette with repeated colors
    my_colors <- rep(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5"), 2)
    
    top_var <- head(summary(pca_results)$importance[2, ], N)
    top_var_percent <- 100 * top_var
    # Create a vector of the top N principal components with their percentage contribution
    pcs <- paste0("PC", seq_along(top_var_percent), " (", round(top_var_percent, 2), "%)")
    
    beeswarm_plot <- biplot_select %>%
      pivot_longer(cols = PC1:N, names_to = "PC", values_to = "value") %>%
      ggplot(aes(x = factor(PC, levels = paste0("PC", 1:N)), y = value, color = PC)) +
      geom_quasirandom(size = 0.6, width = .3) +
      scale_color_manual(values = my_colors, labels = pcs) +
      theme_classic() +
      labs(x = "PCs", y = "Values", title = "Distribution of Values across Multiple PCs ")
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
  
  
  #Generate barplot for top pathways of fgsea
  fgsea_top_pathways <- function(fgsea_results, threshold){
    top_positive_nes <- arrange(fgsea_results, desc(NES)) %>%
      dplyr::filter(padj < 10^(threshold) & NES > 0)
    top_negative_nes <- arrange(fgsea_results, desc(NES)) %>%
      dplyr::filter(padj < 10^(threshold) & NES < 0)
    NES_barplot <- dplyr::bind_rows(top_positive_nes, top_negative_nes)%>%
      ggplot() +
      geom_col(aes(x=reorder(pathway,+NES), y=NES, fill = NES > 0))+
      scale_fill_manual(values =c('TRUE' = 'red', 'FALSE' = 'blue')) +
      theme_minimal() +
      coord_flip()+
      theme(legend.position = "none", axis.text.y = element_text( hjust =1 ,size= 7),axis.title.x = element_text(size = 10))+ 
      labs(title="fgsea results for C2 curated gene sets", x= "",y= "Normalized Enrichment Score (NES)")
    return(NES_barplot)
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
  output$sample_table <- DT::renderDataTable({
    dataf <- sample_data()
    if(is.null(dataf))
      return(null)
    dataf
  }, options = list(
    columnDefs = list(
      list(className = "dt-center", targets = "_all"),
      list(className = "dt-right", targets = "_all", render = JS("function(data, type, full, meta) {return '<div style=\\'width:50px; overflow:hidden\\'>' + data + '</div>';}"))
    )
  ))
  
  #output sample histogram plot
  output$sample_plot <- renderPlot({
    dataf <- sample_data() 
    if(is.null(dataf))
      return(null)
    hist_plot <-histogram(dataf, input$sample_x)
    hist_plot
  })
  
  #output count filter table
  output$filter_count <- renderTable({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    result_tab <- filter_table(dataf, input$slid_var ,input$slid_zero)
    result_tab
  }) 
  
  #output count scatter plot
  output$count_scatter <- renderPlot({
    dataf <- count_data()
    if(is.null(dataf))
      return(null)
    plot1 <- plot_variance_vs_median(dataf,input$slid_var)
    plot2 <- plot_nonzero_vs_median (dataf,input$slid_zero)
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
    beeswarm_plot <- plot_beeswarm(dataf,input$top_PC)
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
  
  #Output download fgsea table
  output$download_fgsea_table <- downloadHandler(
    filename = function() {
      paste("fgsea_table", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(gsea_table(), file, row.names = FALSE)
    }
  )

  observeEvent(input$path_slid, {
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
  })
  
  observeEvent(input$all_path, {
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
  })
  
  #Output fgsea filter scatter plot
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

# Run the application
shinyApp(ui = ui, server = server)
