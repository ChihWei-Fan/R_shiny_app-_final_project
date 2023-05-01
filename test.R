plot_variance_vs_median <- function(data, pass_filter, scale_y_axis=FALSE) {
  new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
    add_column(variance = apply(data[,-1], 1, var), rank = rank(.)) %>%
    mutate(determine = ifelse(variance >= pass_filter, "pass", "fail")) %>%
    ggplot(aes(x= rank, y= variance, color=determine)) +
    geom_point() +
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
    geom_smooth() +
    labs(title= "Median Count vs Variance", x="Rank(Median)", y = "Variance", color="Filter") +
    scale_y_log10()
  return(new_tib)
}

data <- read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.csv", sep="\t", header = TRUE, stringsAsFactors = FALSE)%>%
  as_tibble()%>%
  rename(gene = X)

pass_filter <- 10 # set the threshold for filtering here
pass_filter2 <-5
plot_variance_vs_median(data, pass_filter)
plot_variance_vs_nonzero(data, pass_filter2)


new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
  add_column(variance = apply(data[,-1], 1, var), rank = rank(.)) %>%
  mutate(pass_filter = ifelse(variance >= 10, "pass", "fail"))
new_tib <- tibble(count_median = apply(data[,-1], 1, median)) %>%
  add_column(num_zeros = apply(data[, -1], 1, function(x) sum(x == 0)), rank = rank(.)) %>%
  mutate(determine = ifelse(num_zeros >= pass_filter2, "pass", "fail"))

