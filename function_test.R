library(tidyverse)
library(GEOquery)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(beeswarm)

sample_data <-read.table("GSE64403_FPKM_table.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%as_tibble()   
unknowndata <-read.csv("GSE64810_series_matrix.txt", sep="\t", header = FALSE, stringsAsFactors = FALSE)%>%as_tibble()       

gse <-getGEO(GEO = "GSE64810", GSEMatrix = TRUE)
#gse <-gse[1]
metadata <- pData(phenoData(gse[[1]]))

unknowndata_info <-unknowndata



unknowndata_info <- unknowndata %>% 
  mutate(V1 = apply(unknowndata["V1"],1, function(x) gsub("!", "", x)))%>%
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
names(unknowndata_info) <- as.character(unknowndata_info[1,])
unknowndata_info <- unknowndata_info[-1,] %>%
  mutate(across(c(4,8:11), as.double))

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
  as_tibble()%>%
  rename(gene = X)

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


#
pca_results <- prcomp(scale(t(data[,-1]%>%as.data.frame())), center=FALSE, scale=FALSE)
plot_tibble <- as_tibble(pca_results$x)%>%
  add_column( sample = rownames(pca_results$x), .after=0)
meta <- tibble(sample = rownames(pca_results$x)) %>%
  mutate(Diagnosis = if_else(row_number() <= 49, "normal", "Huntington's Disease"))

biplot <- dplyr::left_join(plot_tibble, meta, by= "sample")
biplot$Diagnosis <- factor(ifelse(biplot$Diagnosis == "normal" & seq_len(nrow(biplot)) <= 49, "normal", "Huntington's Disease"))
biplot <- dplyr::select(biplot,sample, PC1, PC2, Diagnosis)
beeswarm_plot <- ggplot(biplot, aes(x = PC1, y = PC2, color = Diagnosis)) +
  beeswarm(data = biplot, x = biplot[PC1], y = biplot[PC2], pch = 21, col = biplot$Diagnosis) +
  scale_color_manual(values = c("darkorange", "deepskyblue3")) +
  labs(x = "PC1", y = "PC2", color = "Diagnosis") +
  theme_classic()



summary(pca_results)





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


set.seed(123)
group <- sample(c("A", "B", "C"), 100, replace = TRUE)
value <- rnorm(100)
df<-data.frame(group, value)

ggplot(data = data.frame(group, value), aes(x = group, y = value)) +
  geom_quasirandom(alpha = 0.8, color = "black", size = 1) +
  theme_classic() +
  labs(x = "Group", y = "Value")



