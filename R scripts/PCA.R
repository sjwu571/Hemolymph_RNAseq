#### PCA ###

### load filtered combined count matrix #### 
combined_counts_minCPM0.5_nlib1 <- readRDS("./combined_counts_minCPM0.5_nlib1.rds")

### read in metadata ####
sample_names <- read.table("./biosamples_combined.txt", sep = "\t", header = T)
head(sample_names)

sample_table <- data.frame(sample_name = sample_names$sample, 
                           sex = factor(sample_names$sex),
                           time = factor(sample_names$time))
sample_table

### recreate DEseq2 object from count matrix ####
library("DESeq2")
library(ggplot2)
dds <- DESeqDataSetFromMatrix(countData = combined_counts_minCPM0.5_nlib1,
                              colData = sample_table,
                              design = ~ sex)
dds

#### vst, estimate dispersion trend and apply a variance stabilizing transformation####
vsd <- vst(dds)

#### Figure  #####
plotPCA(vsd, "sex") # PC1 vs PC2 
plotPCA(vsd, "sex", pcsToUse = c(1,3))
plotPCA(vsd, "sex", pcsToUse = c(1,4))
plotPCA(vsd, "sex", pcsToUse = c(2,3)) 
# get percent variances

#### making figure ###
library(patchwork)  # For arranging plots

# Common theme settings for all plots
my_theme <- theme_minimal() +
  theme(
    axis.title.x = element_text(size = 9),  # Smaller x-axis label
    axis.title.y = element_text(size = 9),  # Smaller y-axis label
    axis.text = element_text(size = 8),     # Smaller axis tick labels
    legend.title = element_text(size = 9),  # Legend title size
    legend.text = element_text(size = 8),   # Legend text size
    legend.position = "none"                # Remove legend from individual plots
  )

# Panel A: PC1 vs PC2
pcaData <- plotPCA(vsd, intgroup = "sex", pcsToUse = c(1,2), returnData=TRUE) 
panel_a <- ggplot(pcaData, aes(x=PC1, y=PC2, color=sex)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("M"="blue", "F"="orange")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC2: 14% variance") # get these variance values from above

# Panel B: PC1 vs PC3
pcaData <- plotPCA(vsd, intgroup = "sex", pcsToUse = c(1,3), returnData=TRUE) 
panel_b <- ggplot(pcaData, aes(x=PC1, y=PC3, color=sex)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("M"="blue", "F"="orange")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC3: 9% variance")

# Panel C: PC2 vs PC3
pcaData <- plotPCA(vsd, intgroup = "sex", pcsToUse = c(2,3), returnData=TRUE) 
panel_c <- ggplot(pcaData, aes(x=PC2, y=PC3, color=sex)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("M"="blue", "F"="orange")) +
  my_theme +
  labs(x="PC2: 14% variance", y="PC3: 9% variance")

# Panel D: PC1 vs PC4
pcaData <- plotPCA(vsd, intgroup = "sex", pcsToUse = c(1,4), returnData=TRUE) 
panel_d <- ggplot(pcaData, aes(x=PC1, y=PC4, color=sex)) +
  geom_point(size=2.5) +
  scale_color_manual(values=c("M"="blue", "F"="orange")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC4: 8% variance")

# Combine plots with a single centered legend at the bottom
final_figure <- (panel_a + panel_b) / (panel_c + panel_d) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(", tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) +
  plot_layout(guides = "collect") &  # Collect all legends into one
  theme(legend.position = "bottom",  # Place legend at the bottom
        legend.justification = "center")  # Center the legend

print(final_figure)

ggsave(
  "PCA_All_v2.png",
  plot = final_figure,
  width = 6,       # Adjust width (inches)
  height = 6,       # Adjust height
  dpi = 600,        # High resolution for journals
  bg = "white"      # Background color
)




#### Loadings on PC4 #######
mat <- assay(vsd)
# transform the matrix, run PCA
pca <- prcomp(t(mat), scale. = TRUE)

loadings <- pca$rotation[,"PC4"]  # Extract loadings for PC 

sorted_loadings <- sort(loadings, decreasing = TRUE)  # Sort by contribution 

top_genes <- head(sorted_loadings, n = 20)  # Top 20 contributing genes 
top_genes
data.frame(top_genes)


##### color by time (not in the manuscript) #####
# Panel A: PC1 vs PC2
pcaData <- plotPCA(vsd, intgroup = "time", pcsToUse = c(1,2), returnData=TRUE) 
panel_a <- ggplot(pcaData, aes(x=PC1, y=PC2, color=time)) +
  geom_point(size=2.5) +
  scale_color_manual(values = c("T0"="black", "T1"="blue", "T2"="orange", "T3"="red")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC2: 14% variance")

# Panel B: PC1 vs PC3
pcaData <- plotPCA(vsd, intgroup = "time", pcsToUse = c(1,3), returnData=TRUE) 
panel_b <- ggplot(pcaData, aes(x=PC1, y=PC3, color=time)) +
  geom_point(size=2.5) +
  scale_color_manual(values = c("T0"="black", "T1"="blue", "T2"="orange", "T3"="red")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC3: 9% variance")

# Panel C: PC2 vs PC3
pcaData <- plotPCA(vsd, intgroup = "time", pcsToUse = c(2,3), returnData=TRUE) 
panel_c <- ggplot(pcaData, aes(x=PC2, y=PC3, color=time)) +
  geom_point(size=2.5) +
  scale_color_manual(values = c("T0"="black", "T1"="blue", "T2"="orange", "T3"="red")) +
  my_theme +
  labs(x="PC2: 14% variance", y="PC3: 9% variance")

# Panel D: PC1 vs PC4
pcaData <- plotPCA(vsd, intgroup = "time", pcsToUse = c(1,4), returnData=TRUE) 
panel_d <- ggplot(pcaData, aes(x=PC1, y=PC4, color=time)) +
  geom_point(size=2.5) +
  scale_color_manual(values = c("T0"="black", "T1"="blue", "T2"="orange", "T3"="red")) +
  my_theme +
  labs(x="PC1: 23% variance", y="PC4: 8% variance")

# Combine plots with a single centered legend at the bottom
final_figure <- (panel_a + panel_b) / (panel_c + panel_d) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(", tag_suffix = ")",
    theme = theme(
      plot.tag = element_text(size = 10, face = "bold"),
      plot.margin = margin(5, 5, 5, 5)
    )
  ) +
  plot_layout(guides = "collect") &  # Collect all legends into one
  theme(legend.position = "bottom",  # Place legend at the bottom
        legend.justification = "center")  # Center the legend

print(final_figure)
ggsave(
  "PCA_by_Time.png",
  plot = final_figure,
  width = 6,       # Adjust width (inches)
  height = 6,       # Adjust height
  dpi = 600,        # High resolution for journals
  bg = "white"      # Background color
)



