### Volcano plots, between sex and between time points ####

sample_names <- read.table("./biosample_result_A.txt", sep = "\t", header = T)
x <- readRDS("./counts_A_minCPM0.5_nlib1.rds")
dim(x)

# DEG
# Run DESeq2--------------------
library(DESeq2)
FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering, default

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_names,
  design = ~sex
)
dds <- DESeq2::DESeq(dds)

# Extract results--------------------
res <- DESeq2::results(dds, 
                       contrast = c("sex", "F", "M"),
                       independentFiltering = TRUE,
                       alpha = alpha)

# Create volcano plot--------------------
library(ggplot2)
library(ggrepel)

# Convert to data frame and remove NA values
res_df <- as.data.frame(res)
res_df <- res_df[complete.cases(res_df), ]

# Add significance column
res_df$significant <- ifelse(res_df$padj < FDR & abs(res_df$log2FoldChange) > log2(FC), 
                             "Significant", "Not significant")

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not significant" = "gray", "Significant" = "red")) +
  geom_hline(yintercept = -log10(FDR), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-log2(FC)), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = log2(FC), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot: Female vs Male",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)",
       color = "Significance") +
  theme_minimal() +
  theme(legend.position = "top")
             
             # Optional: Add labels for top significant genes
             top_genes <- res_df[res_df$significant == "Significant", ]
             top_genes <- top_genes[order(top_genes$padj), ][1:10, ] # Top 10 most significant
             
             volcano_plot + 
               geom_text_repel(data = top_genes,
                               aes(label = rownames(top_genes)),
                               size = 3,
                               max.overlaps = 20,
                               box.padding = 0.5)
             
#######
# T1
sample_names <- read.table("./biosample_result_B.txt", sep = "\t", header = T)
x <- readRDS("./counts_B_minCPM0.5_nlib1.rds") # then run from line 8 above to make plot, same below
  
# T2           

sample_names <- read.table("./biosample_result_C.txt", sep = "\t", header = T)
x <- readRDS("./counts_C_minCPM0.5_nlib1.rds")

# T3
 
sample_names <- read.table("./biosample_result_D.txt", sep = "\t", header = T)
x <- readRDS("./counts_D_minCPM0.5_nlib1.rds")




# between time
sample_names <- read.table("./biosamples_combined_2.txt", sep = "\t", header = T)
head(sample_names)

# only F or M
sample_table <- subset(sample_names, sex=="M")
dim(sample_table) # 18 bc 2 missing females at T0


#x <- readRDS("./counts_combined_F_only_minCPM0.5_nlib1.rds")
x <- readRDS("./counts_combined_M_only_minCPM0.5_nlib1.rds")

library(DESeq2) # D.E.G.
FC <- 2 # Fold-change cutoff
FDR <- 0.05 # FDR cutoff
alpha <- 0.1 # independent filtering, default

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = x,
  colData = sample_table,
  design = ~ ind + time)

dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, 
                       contrast = c("time", "T3", "T2"), # change as needed
                       independentFiltering = TRUE,
                       alpha = alpha)         

# go to line 28 to draw plot
             