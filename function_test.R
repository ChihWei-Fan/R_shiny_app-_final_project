library(tidyverse)
library(GEOquery)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggbeeswarm)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(DT)
library(biomaRt)
library(fgsea)

sample_data <-read.table("GSE64403_FPKM_table.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%as_tibble()   
       

gse <-getGEO(GEO = "GSE64810", GSEMatrix = TRUE)
#gse <-gse[1]
metadata <- pData(phenoData(gse[[1]]))

unknowndata_info <-unknowndata

read_sample <- function (inputfile){
data <-read.csv(inputfile, sep="\t", header = FALSE, stringsAsFactors = FALSE)%>%as_tibble()
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
}
read<-read_sample("GSE64810_series_matrix.csv")

# histogram

histogram <- function(df, col_name) {
  df[,11]<- log2(df[,11])
  ggplot(df, aes(x = !!sym(col_name))) +
    geom_histogram(binwidth = 1, color = "black", fill = "purple") +
    labs(x = col_name, y = "Count")
}
histogram(read,"Age_of_death")



sidebarLayout(
  sidebarPanel(
    radioButtons(inputId = "x_axis", label = "select y variable", choices = c("Post_mortem_interval","Age_of_death","RNA_integrity_number","mRNA_seq_reads"), selected = "Age_of_death"),
    # Add color inputs
    colourInput(inputId = "base", label = "Choose color", value = "#07B377"),
    #Add a submit buttom
    submitButton(text = "plot",icon = icon("chart-line"))
  ),
  mainPanel(
    plotOutput("sample_plot",width = "90%", height = "550px")
  )
)











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
    else { col_info$mean[i] <- NA #paste(unique(data[[i]]), collapse = "/")
    }
  }
  return(col_info%>%as_tibble())
}
summary_tablef(unknowndata_info) 

#############
data <-read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
  as_tibble()%>%rename(data, gene = X)
colnames(data)[1] <-"gene"

filter_res <- function(data, pass_filter1, pass_filter2) {
  filt_res <- data %>% 
              mutate(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0)),
                     variance = apply(data[,-1], 1, var)) %>%
            filter(variance > 10 & num_zeros > 3)%>% as.data.frame()
  rownames(filt_res) <-filt_res[,1]
  return (filt_res)
  }

filtdata <-filter_res(data,10,3)


plot_heatmap <- function(filter_data) {
  coul <- rev(brewer.pal(11, 'RdBu'))
  num_matrix <- filter_data[-c(1, 71, 72)] %>% as.matrix()%>% log2()
  num_matrix[!is.finite(num_matrix)] <- NA
  heatmap.2(num_matrix, col = coul,trace = "none")
}
plot_heatmap(filtdata)


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
  geom_quasirandom(size = 0.4, width = .3) +
  scale_color_manual(values = my_colors) +
  theme_classic()+
  labs(x = "PCs", y = "Values")
beeswarm_plot
}
plot_beeswarm(data)

summary(pca_results)

##DESeq2
make_se <- function(counts_csv, metafile_csv, subset){
  select_matrix <- read.table(counts_csv, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
    as_tibble()
    rename(gene = X)%>%
    mutate_if(is.integer, as.numeric)%>%
    dplyr::select(gene,contains(subset))
  rowData <- select_matrix["gene"]
  # remove probeset IDs from tibble and turn into a R dataframe so that we can assign row names
  # since tibbles don't support row names
  select_matrix <-as.matrix(dplyr::select(select_matrix, -gene))
  rownames(select_matrix) <- rowData$gene
  
  meta_data <- read.csv(metafile_csv)%>%
    dplyr::select(samplename, timepoint)%>%
    dplyr::filter(timepoint %in% subset)
  meta_data$timepoint <- as.factor(meta_data$timepoint)
  rownames(meta_data) <- colnames(select_matrix)
  # use relevel to set condition
  
  se <- SummarizedExperiment(assays = list(counts=select_matrix)  ,colData = meta_data)
  return(se)
}

#DESeq2
select_matrix <- read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
  as_tibble()
colnames(select_matrix)[1] <-"gene"
select_matrix<-select_matrix%>% mutate_if(is.integer, as.numeric)
rowData <- select_matrix["gene"]
# remove probeset IDs from tibble and turn into a R dataframe so that we can assign row names
# since tibbles don't support row names
select_matrix <- as.matrix(dplyr::select(select_matrix, -gene))
rownames(select_matrix) <- rowData$gene

meta <- data.frame(sample = colnames(data[-1])) %>%
  mutate(Diagnosis = if_else(row_number() <= 49, "normal", "Huntington's Disease"))
meta$Diagnosis <- as.factor(meta$Diagnosis)
rownames(meta) <- colnames(select_matrix)
se <- SummarizedExperiment(assays = list(counts=select_matrix)  ,colData = meta)

#Run  DESeq2  
  dds <- DESeqDataSet(se, design = ~Diagnosis) 
  dds$Diagnosis <- relevel(dds$Diagnosis, ref = "normal")
  dds <- DESeq(dds)
  
  # resultsNames(dds)  #get the names of the result, so we can specify which result we want later
  res <- as.data.frame(results(dds, name="timepoint_vAd_vs_vP0"))
  deseq_res <- list(dds,res)

#Read in DeSeq2 result csv
  
DESeq_res <- read.table("GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
  as_tibble()%>% dplyr::rename(gene= X)
  
#Label Up or down in result data for running gsea
label_res <- function(deseq2_res, padj_threshold){
    # base on p-adj add a new col to specify each row's condition # case_when https://www.statology.org/dplyr-case_when/
    p0_vs_Ad_de<-dplyr::mutate(deseq2_res, volc_plot_status = dplyr::case_when(log2FoldChange > 0 & padj < padj_threshold ~"UP",
                                                      log2FoldChange < 0 & padj < padj_threshold ~"DOWN",
                                                      TRUE ~ "NS"))
  #dplyr::mutate(volc_plot_status =  ifelse(log2FoldChange > 0 & padj < padj_threshold ,
  #                       log2FoldChange < 0 & padj < padj_threshold ~"DOWN") 
  return(p0_vs_Ad_de)
}
labeled_res<-label_res(DESeq_res,.10)
  
#Running GSEA
#process labeled_res
idchange_res <-labeled_res %>%
  dplyr::mutate(gene_id = sub('\\.[0-9]*$', '', labeled_res$gene))%>%   # or you can use seperate(x,sep="\\.", remove....)
  #remove the original gene column and rename the new gene_id column and move it to the first
  dplyr::select(-gene)%>%
  dplyr::rename("gene"="gene_id")%>%
  dplyr::arrange(pvalue)

# subset the gene col in idchange_res into a vector, this will use to convert from id to symbol (as "values")
rank_id_vec <-dplyr::select(idchange_res,gene)%>%
  dplyr::pull()

id_to_hgnc <- function(convert_vector) {
  mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  df<-dplyr::as_tibble(getBM(attributes = c('ensembl_gene_id', "hgnc_symbol"),
                             values     = convert_vector, 
                             mart       = mart))
  return(df)
}
symbol_tib<-id_to_hgnc(rank_id_vec)%>%as.data.frame()

#subset only the hgnc_symbol and filter out the row is empty and duplicate. Turning it to a vector
human_syms <- symbol_tib[,2]%>%
  # recode empty strings "" by NAs
  na_if("") %>%
  # remove NAs
  na.omit() %>%
  # get all the rows with unique symbol
  unique()

#Create the tibble that filter out the symbol that is not in human_syms
filter_tib <- symbol_tib[symbol_tib$hgnc_symbol %in% human_syms, ]
#Combine the tibble idchange_res (with all the values) with the filter_tib(with wanted symbol) using "inner_join() "
combine_res <- inner_join(idchange_res, filter_tib, by=c("gene"="ensembl_gene_id"))

combine_res2 <- dplyr::select(DESeq_res, symbol, log2FoldChange) %>% 
  na.omit() %>%  #make sure no NA in rows
  dplyr::distinct() %>%  # make sure each row is unique
  dplyr::arrange(desc(log2FoldChange))%>%
  tibble::deframe() #deframe() converts two-column data frames to a named vector or list, using the first column as name and the second column as value.

#Read in gmt file
hallmark_pathways_fgsea <- fgsea::gmtPathways("c2.all.v2023.1.Hs.symbols.gmt")  #"c2.cp.v7.5.1.symbols.gmt")
## match the genes provided in the C2 Canonical Pathways gene sets
fgsea_results <-fgsea(hallmark_pathways_fgsea, combine_res2, minSize = 15, maxSize= 500)
#turn fgsea_result into tibble, to let us perform basic exploration easily
fgsea_results <- fgsea_results%>% as.data.frame()
#remove the underline in pathway name
fgsea_results$pathway <-gsub("_"," ",fgsea_results$pathway)
#descending log2foldchange as a ranking metric
fgsea_result_df <-fgsea_results %>% arrange(padj)%>%as_tibble()
# Convert the last column from list to character
fgsea_result_df$leadingEdge <- sapply(fgsea_result_df$leadingEdge, paste, collapse = ",")

# Write the data to a CSV file
write.csv(fgsea_result_df, "fgsea_result.csv", row.names = FALSE)

# generate top pathways bar plot
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
top_pathways (fgsea_result_df,-20)
#read in fgsea csv file
fgseadata <-read.csv("fgsea_result.csv", header = TRUE, stringsAsFactors = FALSE)

#generate fgsea padj filtered pathway table
#Generate padj filtered table in fGSEA tab
gsea_table <- function(dataf, slider, option) {
  filtered_df <- dplyr::filter(dataf, padj < 10^(slider))
  
  if (option == "All") {
    filtered_df <- filtered_df
  } else if (option == "Positive") {
    filtered_df <- dplyr::filter(filtered_df, NES > 0)
  } else if (option == "Negative") {
    filtered_df <- dplyr::filter(filtered_df, NES < 0)
  }
  filtered_df <- filtered_df %>%
    dplyr::mutate(pval = formatC(ifelse(is.na(pval), 0, pval), format = "e"),
                  padj = formatC(ifelse(is.na(padj), 0, padj), format = "e"))
  return(filtered_df)
}
gsea_table(fgseadata,-20,option = "Positive")

filtered_df <- dplyr::filter(fgseadata, padj < 10^(-20))
filtered_df <- dplyr::case_when(
  option == "All" ~ filtered_df,
  option == "Positive" ~ dplyr::filter(filtered_df, NES > 0),
  option == "Negative" ~ dplyr::filter(filtered_df, NES < 0)
) 

fgseadata <-read.csv("fgsea_result.csv", header = TRUE, stringsAsFactors = FALSE)

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
scatter_plot(fgseadata, -20)



####################################################################################################
rownames(unknowndata_info) <- c("GEO accession", "Channel count","Organism_source","Tissue_source","Diagnosis", "Post_mortem_interval","Age_of_death","RNA_integrity_number","mRNA_seq_reads")%>%
  select(unknowndata_info[,-1])

  rename(Sample_names = title, Organism_source = organism_ch1, Tissue_source = characteristics_ch1, Diagnosis = characteristics_ch1.1,
         Post_mortem_interval = characteristics_ch1.2, Age_of_death = characteristics_ch1.3 , RNA_integrity_number = characteristics_ch1.4,
         mRNA_seq_read = characteristics_ch1.5) %>% 
  mutate(Tissue_source = gsub("tissue: ","", Tissue_source), 
         Diagnosis = gsub("diagnosis: ", "",Diagnosis), 
         Post_mortem_interval = gsub("pmi: ", "",Post_mortem_interval),
         Age_of_death = gsub("age of death: ","", Age_of_death), 
         RNA_integrity_number = gsub("rin: ","", RNA_integrity_number), 
         mRNA_seq_read = gsub("mrna-seq reads: ","",mRNA_seq_read),
         across(c(3,7,8,9,10), as.double))   

##
read_data <- function(filename){
  data <-read.table(filename, sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
    mutate_if(is.integer, as.numeric)%>%
    as_tibble()
  return(data)
}

data<- read_data("GSE64810_mlhd_DESeq2_norm_counts_adjust.csv")%>%
              rename(gene = X)
data <-read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%as_tibble()   

##
filter_genes <- function(counts, var_percentile, nonzero_samples) {
  # Calculate the variance across samples for each gene
  gene_variances <- apply(counts[,-1], 1, var)
  
  # Calculate the percentile of variance above which genes will be kept
  var_threshold <- quantile(gene_variances, probs = var_percentile)
  
  # Find the genes with variance above the threshold
  var_filtered_genes <- gene_variances >= var_threshold
  
  # Count the number of nonzero samples for each gene
  nonzero_counts <- apply(counts[,-1], 1, function(x) sum(x != 0))
  
  # Find the genes with at least nonzero_samples non-zero samples
  nonzero_filtered_genes <- nonzero_counts >= nonzero_samples
  
  # Combine the two sets of filtered genes
  filtered_genes <- var_filtered_genes & nonzero_filtered_genes
  
  # Use the filtered gene list to subset the original matrix
  filtered_matrix <- as_tibble(counts[filtered_genes,])
  
  return(filtered_matrix)
}


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
filter_table(data,5)


    data <- read.csv("GSE64810_series_matrix.csv", sep="\t", header=FALSE, stringsAsFactors=FALSE) %>%
      as_tibble()
    data <- data%>% mutate(V1 = apply(data["V1"],1, function(x) gsub("!", "", x)))%>%
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
    
    run_gsea <- function(labeled_results, gmt, min_size, max_size) {
      
      ##Create a tibble with new gene id (without version) and, other cols of stat values and volc_status col.
      idchange_res <-labeled_results %>%
        dplyr::mutate(gene_id = sub('\\.[0-9]*$', '', labeled_results$gene))%>%   # or you can use seperate(x,sep="\\.", remove....)
        #remove the original gene column and rename the new gene_id column and move it to the first
        dplyr::select(-gene)%>%
        dplyr::rename("gene"="gene_id")%>%
        dplyr::arrange(pvalue)
      
      # subset the gene col in idchange_res, this will use to convert from id to symbol
      rank_id_vec <-dplyr::select(idchange_res,gene)%>%
        dplyr::pull() #turn the ranked ensembl id  tibble into a simple character vector
      
      ## Connect to the biomart                                                               
      human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl',host = "https://dec2021.archive.ensembl.org/")
      
      ##convert the mouse ensembl gene IDs to human HGNC symbols and turn the result to a tibble. use getLDS() # 
      
      genesV2 <- getLDS(mart = mouse, attributes = c('mgi_symbol','ensembl_gene_id'), ##Attributes define the data we are interested in retrieving.(here we want to retrieve the gene symbols and their snsembl_id)  
                        martL = human,attributesL = c('hgnc_symbol','ensembl_gene_id'),
                        #filters = 'ensembl_gene_id', ##Filters and values are used to define restrictions on the query ##https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#how-to-build-a-biomart-query
                        values = rank_id_vec) %>% as_tibble()
      #subset only the hgnc_symbol and filter out the row is empty and duplicate. Turning it to a vector
      human_syms <- genesV2[, 3] %>%
        # recode empty strings "" by NAs
        na_if("") %>%
        # remove NAs
        na.omit() %>%
        # get all the rows with unique symbol
        distinct()
      human_syms <- as.vector(human_syms$HGNC.symbol )
      
      #Create the tibble that filter out the symbol that is not in human_syms
      filter_tib <- genesV2[genesV2$HGNC.symbol %in% human_syms, ]%>%
        dplyr::select(Gene.stable.ID,HGNC.symbol)
      
      #Combine the tibble idchange_res (with all the values) with the filter_tib(with wanted symbol) using "inner_join() "
      combine_res <- inner_join(idchange_res, filter_tib, by=c("gene"="Gene.stable.ID"))
      # only need gene symbol and the log2FC 
      combine_res2 <- dplyr::select(combine_res,HGNC.symbol, log2FoldChange) %>% 
        na.omit() %>%  #make sure no NA in rows
        dplyr::distinct() %>%  # make sure each row is unique
        dplyr::arrange(desc(log2FoldChange))%>%
        tibble::deframe() #deframe() converts two-column data frames to a named vector or list, using the first column as name and the second column as value.
      
      ## read in gene set "c2.cp.v7.5.1.symbols.gmt"
      hallmark_pathways_fgsea <- fgsea::gmtPathways(gmt)  #"c2.cp.v7.5.1.symbols.gmt")
      
      ## match the genes provided in the C2 Canonical Pathways gene sets
      fgsea_results <-fgsea(hallmark_pathways_fgsea, combine_res2, minSize = min_size, maxSize= max_size)
      
      #turn fgsea_result into tibble, to let us perform basic exploration easily
      fgsea_results <- fgsea_results%>% as_tibble()
      #descending log2foldchange as a ranking metric
      fgsea_result_tibble <-fgsea_results%>% arrange(NES)
      
      return(fgsea_result_tibble)
    }
    
    


    # Add an observer to monitor the changes made to input$path_slid
    observeEvent(input$path_slid, {
      # Call gsea_table() with the new input$path_slid value and update the output
      dataf <- fgsea_data()
      filter_tab <- gsea_table(dataf, input$path_slid, input$all_path)
      output$fgsea_filt_table <- DT::renderDataTable({
        filter_tab
      })
    })
    
    
##old read sample
    data <- read.csv(input$sample_file$datapath, sep="\t", header=FALSE, stringsAsFactors=FALSE) %>%
      as_tibble()
    data <- data %>% mutate(V1 = apply(data["V1"],1, function(x) gsub("!", "", x)))%>% as.data.frame()
    slice(c(1:3,6,8:14)) %>%
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

