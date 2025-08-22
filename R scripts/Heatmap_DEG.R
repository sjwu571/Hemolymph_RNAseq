
#### Heatmap DEG #####

library(DESeq2)
library(pheatmap)
library(viridis)
library(ggplot2)

#### T0 ####

sample_names <- read.table("./biosample_result_A.txt", sep = "\t", header = T)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
sample_table

library(DESeq2)
x <- readRDS("./counts_A_minCPM0.5_nlib1.rds")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~sex )

vsd <- vst(dds, blind = TRUE) 
# Extract the normalized counts
normalized_counts <- assay(vsd)

# get the DEG's
DEG <- read.table("./DEG_table_fdr05.txt")

# extract gene id for DEGs
selected_genes <- rownames(DEG)

vsd_subset <- normalized_counts[selected_genes, ]
colnames(vsd_subset) <- sample_table$ind
plot <- pheatmap(vsd_subset,
         scale = "row",  # Scale by row (genes) to highlight relative expression
         cluster_rows = FALSE,  # Cluster rows (genes)
         cluster_cols = TRUE,  # Cluster columns (samples)
         color = viridis(100),  # Color scheme
         show_rownames = FALSE, 
         clustering_method = "average",
         main = "T0 DEG")

ggsave("T0_DEG_heatmap_2.pdf", plot, width = 3.5, height = 4, dpi = 600)

# T1 
sample_names <- read.table("./biosample_result_B.txt", sep = "\t", header = T)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
sample_table

library(DESeq2)
x <- readRDS("~/counts_B_minCPM0.5_nlib1.rds")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~sex )

vsd <- vst(dds, blind = TRUE) 
# Extract the normalized counts
normalized_counts <- assay(vsd)

DEG <- read.table("DEG_table_fdr05.txt")

# extract gene id for DEGs
selected_genes <- rownames(DEG)

vsd_subset <- normalized_counts[selected_genes, ]
colnames(vsd_subset) <- sample_table$ind
plot <- pheatmap(vsd_subset,
                 scale = "row",  # Scale by row (genes) to highlight relative expression
                 cluster_rows = FALSE,  # Cluster rows (genes)
                 cluster_cols = TRUE,  # Cluster columns (samples)
                 color = viridis(100),  # Color scheme
                 show_rownames = FALSE, 
                 clustering_method = "average", #"single" looks better for this time point
                 main = "T1 DEG")


ggsave("T1_DEG_heatmap_2.pdf", plot, width = 3.5, height = 4, dpi = 600)


# T2
sample_names <- read.table("./biosample_result_C.txt", sep = "\t", header = T)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
sample_table

library(DESeq2)
x <- readRDS("./counts_C_minCPM0.5_nlib1.rds")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~sex )


vsd <- vst(dds, blind = TRUE) 
# Extract the normalized counts
normalized_counts <- assay(vsd)

DEG <- read.table("DEG_table_fdr05.txt")

# extract gene id for DEGs
selected_genes <- rownames(DEG)

vsd_subset <- normalized_counts[selected_genes, ]
colnames(vsd_subset) <- sample_table$ind
plot <- pheatmap(vsd_subset,
                 scale = "row",  # Scale by row (genes) to highlight relative expression
                 cluster_rows = FALSE,  # Cluster rows (genes)
                 cluster_cols = TRUE,  # Cluster columns (samples)
                 color = viridis(100),  # Color scheme
                 show_rownames = FALSE, 
                 border_color = NA, # this time point is weird 
                 clustering_method = "average",
                 main = "T2 DEG")

# why the grey border between rows?
ggsave("T2_DEG_heatmap_2.pdf", plot, width = 3.5, height = 4, dpi = 600)

# T3

sample_names <- read.table("./biosample_result_D.txt", sep = "\t", header = T)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time),
                           ind = factor(sample_names$ind))
sample_table

library(DESeq2)
x <- readRDS("./counts_D_minCPM0.5_nlib1.rds")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~sex )

vsd <- vst(dds, blind = TRUE) 
# Extract the normalized counts
normalized_counts <- assay(vsd)

DEG <- read.table("DEG_table_fdr05.txt")

# extract gene id for DEGs
selected_genes <- rownames(DEG)

vsd_subset <- normalized_counts[selected_genes, ]
colnames(vsd_subset) <- sample_table$ind
plot <- pheatmap(vsd_subset,
                 scale = "row",  # Scale by row (genes) to highlight relative expression
                 cluster_rows = FALSE,  # Cluster rows (genes)
                 cluster_cols = TRUE,  # Cluster columns (samples)
                 color = viridis(100),  # Color scheme
                 show_rownames = FALSE, 
                 clustering_method = "average",
                 main = "T3 DEG")

# why the grey border between rows?
ggsave("T3_DEG_heatmap_2.pdf", plot, width = 3.5, height = 4, dpi = 600)

