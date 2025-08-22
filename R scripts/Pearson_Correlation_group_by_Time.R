##### correlation #####

combined_counts_minCPM0.5_nlib1 <- readRDS("./combined_counts_minCPM0.5_nlib1.rds")

x <- combined_counts_minCPM0.5_nlib1

sample_names <- read.table("./biosamples_combined.txt", sep = "\t", header = T)
head(sample_names)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time), 
                           ind = factor(sample_names$ind))
sample_table

### Remake dds file to do vst tranformation ###

library(DESeq2)
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~ sex + time) 

# Perform VST normalization
vsd <- vst(dds, blind = TRUE)  # blind = TRUE for unsupervised normalization


# Extract the normalized counts
normalized_counts <- assay(vsd)

# View the first few rows of the normalized counts
head(normalized_counts)


#### calculate correlation ####
correlation_matrix <- cor(normalized_counts, method = "pearson")

# View the correlation matrix
colnames(correlation_matrix)
correlation_matrix
mean(correlation_matrix) #0.9286836

# Install pheatmap if not already installed
# if (!requireNamespace("pheatmap", quietly = TRUE)) {
#   install.packages("pheatmap")
# }
library(pheatmap)
library(dplyr)
library(ggplot2)

#### arrange samples first by sex, then time####
df_sorted <- sample_table %>% arrange (time, sex)
head(df_sorted)
tail(df_sorted)


# Ensure the correlation matrix rows and columns are ordered by sample_name in df_sorted
ordered_samples <- df_sorted$sample_name
correlation_matrix_ordered <- correlation_matrix[ordered_samples, ordered_samples]


#### try different color#####
library(viridis)
pheatmap(correlation_matrix_ordered,
         cluster_rows = FALSE,  # Cluster rows based on correlation
         cluster_cols = FALSE,  # Do not cluster columns
         display_numbers = FALSE,
         number_format = "%.2f",
         color = viridis(100),
         main = "")


# to highlight samples with highest corr efficient #
pheatmap(correlation_matrix_ordered,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         display_numbers = FALSE,
         number_format = "%.2f",
         color = viridis(100),  
         breaks = seq(0.9, 1, length.out = 100), 
         clustering_method = "average", # or default
         main = "")


#### rename colnames and rownames on heat map ####

# View the correlation matrix
colnames(correlation_matrix_ordered) <- df_sorted$ind
rownames(correlation_matrix_ordered) <- df_sorted$ind

# to highlight samples with highest corr efficient #
plot <- pheatmap(correlation_matrix_ordered,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         display_numbers = FALSE,
         number_format = "%.2f",
         color = viridis(100),  
         breaks = seq(0.9, 1, length.out = 100), 
         clustering_method = "average") # or default
      #   main = "Pearson Correlation Matrix (Vst Normalized Data)")

ggsave("Pearson_heatmap_by_time.pdf", plot , width = 8, height = 6, dpi = 600)

